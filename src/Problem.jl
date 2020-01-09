mutable struct Parameters
    n1fctrs::Tuple{Real,Real}
    n2ton1fctrs::Tuple{Real,Real}
end

struct Problem
    grids
    n1fctrs
    n1values
    n2mincontinueabs
    n2ton1fctrs
    maxmultipleonestage
    nmax
    curtail_stage_one_buffer
    x1min
    partial
    unimodal
    type
    objective
    power
    toer

    function Problem(
        objective,
        toer,
        power;
        type                  = :TwoStage,
        unimodal              = false,
        partial::Tuple{TI,TI} = (0, 0),
        maxmultipleonestage   = 2.,
        α                     = toer.α,
        β                     = power.β,
        pnull                 = mean(toer.score.prior),
        palt                  = mean(power.score.prior),
        nmax                  = min(
                150,
                maxmultipleonestage*one_stage_sample_size(pnull, α, palt, β)
            ) |> ceil |> Int,
        n1fctrs::Tuple{Real,Real} = (.2/maxmultipleonestage, .51),
        n1min                 = max(partial[2], n1fctrs[1] * nmax) |> ceil |> Int,
        n1max                 = n1fctrs[2] * nmax |> ceil |> Int,
        n1values              = collect(n1min:n1max),
        n2mincontinueabs      = 5,
        n2ton1fctrs::Tuple{Real,Real} = (1.1, 4.0),
        curtail_stage_one_buffer = 2
    ) where {TI<:Integer}

        @assert 0 <= partial[1] <= partial[2]

        x1min = partial[1]

        # prebuild sparse grid space

        function n2(n1, x1)

            if type == :OneStage; return [0] end
            n2min = max(n2mincontinueabs, n1*(n2ton1fctrs[1] - 1)) |> ceil |> Int
            n2max = min(nmax - n1, n1*(n2ton1fctrs[2] - 1)) |> floor |> Int
            # curtail if error rates are too low in stage one
            if (1 - cdf(x1 - curtail_stage_one_buffer, n1, toer.score.prior; partial = partial) <= α) |
                (cdf(x1 + curtail_stage_one_buffer, n1, power.score.prior; partial = partial) <= β)
                return [0]
            end
            # add 0 for early stopping
            return n2min > 0 ? vcat( [0], collect(n2min:n2max) ) : collect(n2min:n2max)
        end

        function c2(n1, x1, n2)

            if n2 == 0; return [-Inf, Inf] end
            candidates = collect(0:(n2 - 1))
            # filter Pr[X1 == x1] * toer[x1] > alpha
            prior      = toer.score.prior
            candidates = candidates[ pmf(x1, n1, prior; partial = partial) .* (1 .- cdf.(candidates, n2, update(prior, x1, n1))) .<= α ]
            # filter Pr[X1 == x1] * tter[x1]) > beta
            prior      = power.score.prior
            candidates = candidates[ pmf(x1, n1, prior; partial = partial) .* cdf.(candidates, n2, update(prior, x1, n1)) .<= β ]
            # check conditional power and type one error rate
            valid = trues(length(candidates))
            for i in 1:length(candidates)
                for cnstr in (power, toer)
                    min, max = cnstr.conditional
                    valid[i] = !(min <= (1 - cdf(candidates[i], n2, update(cnstr.score.prior, x1, n1))) <= max) ? false : valid[i]
                end
            end
            return vcat([-Inf], candidates[valid], [Inf]) # required to be safe!
        end

        function n1_grid(n1)
            [ (x1, n2, c2)
                for x1 in x1min:n1
                for n2 in n2(n1, x1)
                for c2 in c2(n1, x1, n2)
            ]
        end

        grids = Dict(n1 => n1_grid(n1) for n1 in n1values)

        return new(
            grids,
            n1fctrs,
            sort(n1values),
            n2mincontinueabs,
            n2ton1fctrs,
            maxmultipleonestage,
            nmax,
            curtail_stage_one_buffer,
            x1min,
            partial,
            unimodal,
            type,
            objective,
            power,
            toer
        )
    end

end

Base.iterate(problem::Problem, state = 0) = state > 0 ? nothing : (problem, state + 1)
Base.length(problem::Problem) = 1

Base.show(io::IO, problem::Problem) = print(io, string(problem))
Base.show(io::IO, ::MIME"application/prs.juno.inline", problem::Problem) = print(io, string(problem))

Base.string(problem::Problem) = @sprintf "Problem<x1>=%i,%i<=n1<=%i,n<=%i>" problem.x1min problem.n1values[1] problem.n1values[end] problem.nmax

# query precomputed grid values
grid(problem::Problem, n1)  = problem.grids[n1]
grid(problem::Problem)      = [[ (n1, x1, n2, c2) for (x1, n2, c2) in grid] for (n1, grid) in problem.grids] |> x -> vcat(x...)

Base.size(problem::Problem, n1) = size(problem.grids[n1],1)
Base.size(problem::Problem) = [size(problem, n1) for n1 = problem.n1values] |> sum



function build_model(problem::Problem; verbose::Bool = true)

    n1values = problem.n1values
    if verbose
        prog = ProgressMeter.Progress(
            1 + length(n1values) + 3,
            desc      = "Building ILP Problem: ",
            dt        = 0.1,
            barglyphs = ProgressMeter.BarGlyphs("[=> ]"),
            barlen    = 20,
            color     = :gray
        )
    end
    update_progress(str) = if verbose; ProgressMeter.next!(prog; showvalues = [(Symbol(str), "...")], valuecolor = :gray) end

    update_progress("initialising JuMP variables")
    m = Model()
    # ind[x1, n1, n2, c2] == 1 iff n1 = n1, n2(x1) = n2, c2(x1) = c2
    @variable(m, ind[(n1, x1, n2, c2) in grid(problem)], Bin)
    # need to make sure that exactly one n1 value is selected, IP or trick via
    # auxiliary variables: n1_selected[n1] == 1 iff n1 = n1
    @variable(m, n1_selected[n1 in n1values], Bin)
    @constraint(m, sum(n1_selected[n1] for n1 in n1values) == 1)
    if (problem.type == :TwoStage) & problem.unimodal
        @variable(m, is_mode[n1 in n1values, x1 in collect(problem.x1min:n1)], Bin)
        for n1_ in n1values
            # at least one mode (can be on boundary as well!)
            @constraint(m,
              sum( is_mode[n1_, x1] for x1 in collect(problem.x1min:n1_) ) >= 1
            )
        end
    end
    # implement all constraints conditional on n1
    for n1_ in n1values
        update_progress(@sprintf "building constraints for n1=%i" n1_)
        # make sure that n1_selected[n1] is 1 iff any of the other values is assigned
        @constraint(m,
            sum( ind[(n1_, x1, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) )
            <= 5*problem.nmax*n1_selected[n1_]
        )
        for x1_ in collect(problem.x1min:n1_)
            # make it a function in x1
            @constraint(m,
                sum( ind[(n1_, x1_, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) if x1 == x1_ )
                == n1_selected[n1_]
            )
            # make sure that n1(x1) is constant
            if problem.x1min < x1_
                @constraint(m,
                    sum( ind[(n1_, x1_, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) if (x1 == x1_) )
                    == sum( ind[(n1_, x1_ - 1, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) if (x1 == x1_ - 1) )
                )
                if problem.type == :GroupSequential
                    for n2_ in unique([n2 for (x1, n2, c2) in grid(problem, n1_) if n2 > 0])
                        @constraint(m,
                            sum( ind[(n1_, x1_, n2_, c2)] # same n2
                                for (x1, n2, c2) in grid(problem, n1_)
                                if (x1 == x1_) & (n2 == n2_) & isfinite(c2)
                            ) +
                            ind[(n1_, x1_, 0, -Inf)] # or early efficacy
                            >= sum( ind[(n1_, x1_ - 1, n2_, c2)]
                                for (x1, n2, c2) in grid(problem, n1_)
                                if (x1 == x1_ - 1) & (n2 == n2_) & isfinite(c2)
                            )
                        )
                    end
                end
            end
            # make sure we have contiguous stopping for futility ...
            if x1_ > problem.x1min
                @constraint(m, ind[(n1_, x1_ - 1, 0, Inf)] >= ind[(n1_, x1_, 0, Inf)] )
            end
            # ... and for efficacy
            if x1_ < n1_
                @constraint(m, ind[(n1_, x1_ + 1, 0, -Inf)] >= ind[(n1_, x1_, 0, -Inf)] )
            end
            # optional unimodality
            if (problem.type == :TwoStage) & problem.unimodal
                for xx1_ in (problem.x1min + 1):x1_
                    # for a mode at x1_, n must be increasing before ...
                    @constraint(m,
                           sum( n2 * ind[(n1_, xx1_, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) if (x1 == xx1_) )
                         - sum( n2 * ind[(n1_, xx1_ - 1, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) if (x1 == xx1_ - 1) )
                        >= 10*problem.nmax*(is_mode[n1_, x1_] - 1)
                    )
                end
                for xx1_ in (x1_ + 1):n1_
                    # ... and decreasing after
                    @constraint(m,
                           sum( n2 * ind[(n1_, xx1_, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) if (x1 == xx1_) )
                         - sum( n2 * ind[(n1_, xx1_ - 1, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) if (x1 == xx1_ - 1) )
                        <= -10*problem.nmax*(is_mode[n1_, x1_] - 1)
                    )
                end
            end
        end
    end
    update_progress("adding type one error rate constraint")
    @constraint(m,
        sum( integrand_x1(problem.toer, x1, n1, n2, c2; partial_stage_one = problem.partial) * ind[(n1, x1, n2, c2)]
            for (n1, x1, n2, c2) in grid(problem)
        ) <= problem.toer.α
    )
    update_progress("adding power constraint")
    @constraint(m,
        sum( integrand_x1(problem.power, x1, n1, n2, c2; partial_stage_one = problem.partial) * ind[(n1, x1, n2, c2)]
            for (n1, x1, n2, c2) in grid(problem)
        ) >= 1 - problem.power.β
    )
    update_progress("adding objective")
    add!((m, ind), problem.objective, problem; partial = problem.partial)
    update_progress("ready to start ILP solver ")
    if verbose; ProgressMeter.next!(prog) end
    return m, ind, n1_selected
end


function optimise!(m, verbosity::Integer, timelimit::Integer)

    optimize!(
        m,
        with_optimizer(GLPK.Optimizer, msg_lev = verbosity, tm_lim = 1000*timelimit)
    )
    # termination_status(model.jump_model)
end


function extract_solution(problem::Problem, ind, n1_selected)

    vals = value.(ind)
    # find n1
    x1obs = problem.partial[1]
    n1_ = problem.n1values[findfirst(value.(n1_selected).data .== 1)]
    c2_res = repeat([Inf], n1_ + 1 - x1obs)
    n2_res = zeros(n1_ + 1 - x1obs)
    for (x1, n2, c2) in grid(problem, n1_)
        if vals[(n1_, x1, n2, c2)] == 1
            c2_res[x1 + 1 - x1obs] = c2
            n2_res[x1 + 1 - x1obs] = n2
        end
    end
    return collect(x1obs:n1_), n2_res, c2_res
end



mutable struct OptimalDesign{TI<:Integer,TR<:Real} <: AbstractDesign
    n2::Vector{TI}
    c2::Vector{TR}
    problem::Problem
    info::Dict{String,Any}
end
function OptimalDesign(design::Design, problem::Problem, info::Dict{String,Any})
    return OptimalDesign{eltype(design.n2),eltype(design.c2)}(design.n2, design.c2, problem, info)
end

function optimise(problem::Problem; verbosity = 3, timelimit = 180)

    info = Dict{String,Any}()
    _, info["model build time [s]"], _, _, _ = @timed begin
        m, ind, n1_selected = build_model(problem; verbose = verbosity > 0)
    end
    _, info["model ILP solution time [s]"], _, _, _ = @timed begin
        optimise!(m, verbosity, timelimit)
    end
    _, info["solution extraction time [s]"], _, _, _ = @timed begin
        xx1, nn2, cc2 = extract_solution(problem, ind, n1_selected)
    end
    info["total time [s]"] = sum([val for (key, val) in info])
    info["number of variables"] = size(problem)
    if problem.partial[2] > 0
        return xx1, nn2, cc2, info
    else
        design  = Design(nn2, cc2)
        pnull   = problem.toer.score.bounds[2]
        α       = problem.toer.α
        mlecomp = mlecompatible(design, pnull, α)
        if !mlecomp["compatible"] & (verbosity > 0)
            @warn @sprintf "design is not compatible with MLE-ordering, incompatibility degree is %i/%i" mlecomp["incompatibility degree"] size(mlecomp["details"], 1)
        end
        info["MLE-compatible"] = mlecomp
        return OptimalDesign(design, problem, info)
    end
end



function adapt(design::OptimalDesign, prior::TP, partial::Tuple{TI,TI};
        α         = design.problem.toer.score(design, partial_stage_one = partial),
        β         = 1 - design.problem.power.score(design, partial_stage_one = partial),
        verbosity = 0,
        timelimit = 300
    ) where {TP<:Prior,TI<:Integer}

        objective   = update(design.problem.objective, prior)
        toer        = design.problem.toer
        toer.α      = α
        power       = update(design.problem.power, prior)
        power.β     = β
        pnull, palt = mean(toer.score.prior), mean(power.score.prior)

        n1fctrs     = design.problem.n1fctrs
        n2mincontinueabs = design.problem.n2mincontinueabs
        n2ton1fctrs = design.problem.n2ton1fctrs
        maxmultipleonestage = design.problem.maxmultipleonestage
        curtail_stage_one_buffer = design.problem.curtail_stage_one_buffer
        type = design.problem.type
        unimodal = design.problem.unimodal

        nmax = max(
            min(
                150,
                maxmultipleonestage * one_stage_sample_size(pnull, α, palt, β)
            ) |> ceil |> Int,
            design.problem.nmax
        )
        n1min = max(partial[2], n1fctrs[1] * nmax) |> floor |> Int
        n1max = n1fctrs[2] * nmax |> ceil |> Int

        adaptation_problem = Problem(
            objective,
            toer,
            power;
            partial               = partial,
            maxmultipleonestage   = maxmultipleonestage,
            α                     = α,
            β                     = β,
            pnull                 = pnull,
            palt                  = palt,
            nmax                  = nmax,
            n1fctrs               = n1fctrs,
            n1min                 = n1min,
            n1max                 = n1max,
            n1values              = collect(n1min:n1max),
            n2mincontinueabs      = n2mincontinueabs,
            n2ton1fctrs           = n2ton1fctrs,
            curtail_stage_one_buffer = curtail_stage_one_buffer,
            type                  = design.problem.type,
            unimodal              = unimodal
        )

        m, ind, n1_selected = build_model(adaptation_problem; verbose = verbosity > 0)
        optimise!(m, verbosity, timelimit)
        xx1, nn2, cc2 = extract_solution(adaptation_problem, ind, n1_selected)
        return xx1, nn2, cc2
end
