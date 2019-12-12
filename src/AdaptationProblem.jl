struct AdaptationProblem
    grids
    n1values
    nmax
    n2mincontinueabs
    n2mincontinuereln1
    n2maxcontinuereln1
    x1min
    type
    objective
    power
    toer
    curtail_stage_one_fct
    old_design
    x1partial
    n1partial

    function AdaptationProblem(
            design::OptimalDesign,
            x1partial::TI,
            n1partial::TI
        ) where {TI<:Integer}

        @assert 0 <= x1partial <= n1partial

        n1values              = design.problem.n1values[design.problem.n1values .>= n1partial]
        nmax                  = design.problem.nmax
        n2mincontinueabs      = design.problem.n2mincontinueabs
        n2mincontinuereln1    = design.problem.n2mincontinuereln1
        n2maxcontinuereln1    = design.problem.n2maxcontinuereln1
        x1min                 = x1partial
        type                  = design.problem.type
        objective             = design.problem.objective
        power                 = design.problem.power
        toer                  = design.problem.toer
        curtail_stage_one_fct = design.problem.curtail_stage_one_fct

        # prebuild sparse grid space
        function n2(n1, x1)

            if type == :OneStage; return [0] end
            n2min = max(n2mincontinueabs, Int(ceil(n1*n2mincontinuereln1)) - n1)
            n2max = min(nmax - n1, Int(floor(n1*(n2maxcontinuereln1))) - n1)
            if x1 > 0 # compute stage one toer for rejecting at x1
                toer1 = 1 - cdf(x1 - 1, n1, toer.score.prior; xpartial = x1partial, npartial = n1partial)
                if (toer1 <= curtail_stage_one_fct*toer.α)
                    return [0]
                end
            end
            if x1 < n1 # compute stage one tter for rejecting at x1
                tter1 = cdf(x1 + 1, n1, power.score.prior; xpartial = x1partial, npartial = n1partial)
                if (tter1 <= curtail_stage_one_fct*power.β)
                    return [0]
                end
            end
            # add 0 for early stopping
            res = n2min > 0 ? vcat( [0], collect(n2min:n2max) ) : collect(n2min:n2max)
            @assert res[1] == 0
            return res
        end

        function c2(n1, x1, n2)

            if n2 == 0; return [-Inf, Inf] end
            candidates = collect(0:(n2 - 1))
            # filter Pr[X1 == x1] * toer[x1] > alpha
            prior      = toer.score.prior
            candidates = candidates[ pmf(x1, n1, prior; xpartial = x1partial, npartial = n1partial) .* (1 .- cdf.(candidates, n2, update(prior, x1, n1))) .<= toer.α ]
            # filter Pr[X1 == x1] * tter[x1]) > beta
            prior      = power.score.prior
            candidates = candidates[ pmf(x1, n1, prior; xpartial = x1partial, npartial = n1partial) .* cdf.(candidates, n2, update(prior, x1, n1)) .<= power.β ]
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
            sort(n1values),
            nmax,
            n2mincontinueabs, n2mincontinuereln1, n2maxcontinuereln1,
            x1min,
            type,
            objective,
            power,
            toer,
            curtail_stage_one_fct,
            design,
            x1partial,
            n1partial
        )
    end

end

Base.iterate(problem::AdaptationProblem, state = 0) = state > 0 ? nothing : (problem, state + 1)
Base.length(problem::AdaptationProblem) = 1

Base.show(io::IO, problem::AdaptationProblem) = print(io, string(problem))
Base.show(io::IO, ::MIME"application/prs.juno.inline", problem::AdaptationProblem) = print(io, string(problem))

Base.string(problem::AdaptationProblem) = @sprintf "AdaptationProblem<x1>=%i,%i<=n1<=%i,n<=%i>" problem.x1min problem.n1values[1] problem.n1values[end] problem.nmax

# query precomputed grid values
grid(problem::AdaptationProblem, n1)  = problem.grids[n1]
grid(problem::AdaptationProblem)      = [[ (n1, x1, n2, c2) for (x1, n2, c2) in grid] for (n1, grid) in problem.grids] |> x -> vcat(x...)
x1(problem::AdaptationProblem, n1)    = collect((problem.x1min):n1)

Base.size(problem::AdaptationProblem, n1) = size(problem.grids[n1],1)
Base.size(problem::AdaptationProblem) = [size(problem, n1) for n1 = problem.n1values] |> sum



function build_model(problem::AdaptationProblem; verbose::Bool = true)

    n1values = problem.n1values
    if verbose
        prog = ProgressMeter.Progress(
            1 + length(n1values) + 3,
            desc      = "Building IP AdaptationProblem: ",
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
    # implement all constraints conditional on n1
    for n1_ in n1values
        update_progress(@sprintf "building constraints for n1=%i" n1_)
        # make sure that n1_selected[n1] is 1 iff any of the other values is assigned
        @constraint(m,
            sum( ind[(n1_, x1, n2, c2)] for (x1, n2, c2) in grid(problem, n1_) )
            <= 5*problem.nmax*n1_selected[n1_]
        )
        for x1_ in x1(problem, n1_)
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
        end
        # implement conditional error constraints
        

    end
    update_progress("adding type one error rate constraint")
    @constraint(m,
        sum( integrand_x1(problem.toer, x1, n1, n2, c2) * ind[(n1, x1, n2, c2)]
            for (n1, x1, n2, c2) in grid(problem)
        ) <= problem.toer.α
    )
    update_progress("adding power constraint")
    @constraint(m,
        sum( integrand_x1(problem.power, x1, n1, n2, c2) * ind[(n1, x1, n2, c2)]
            for (n1, x1, n2, c2) in grid(problem)
        ) >= 1 - problem.power.β
    )
    update_progress("adding objective")
    add!((m, ind), problem.objective, problem)
    update_progress("ready to start ILP solver ")
    if verbose; ProgressMeter.next!(prog) end
    return m, ind, n1_selected
end





function extract_solution(problem::AdaptationProblem, ind, n1_selected)

    vals = value.(ind)
    # find n1
    n1_ = problem.n1values[findfirst(value.(n1_selected).data .== 1)]
    c2_res = repeat([Inf], n1_ + 1 - problem.x1partial)
    n2_res = zeros(n1_ + 1 - problem.x1partial)
    for (x1, n2, c2) in grid(problem, n1_)
        if vals[(n1_, x1, n2, c2)] == 1
            c2_res[x1 + 1 - problem.x1partial] = c2
            n2_res[x1 + 1 - problem.x1partial] = n2
        end
    end
    return collect(problem.x1partial:n1_), n2_res, c2_res
end



# mutable struct OptimalDesign{TI<:Integer,TR<:Real} <: AbstractDesign
#     n2::Vector{TI}
#     c2::Vector{TR}
#     problem::AdaptationProblem
#     info::Dict{String,Any}
# end
# function OptimalDesign(design::Design, problem::AdaptationProblem, info::Dict{String,Any})
#     return OptimalDesign{eltype(design.n2),eltype(design.c2)}(design.n2, design.c2, problem, info)
# end
#
# function optimise(problem::AdaptationProblem; verbosity = 3, timelimit = 180)
#
#     info = Dict{String,Any}()
#     _, info["model build time [s]"], _, _, _ = @timed begin
#         m, ind, n1_selected = build_model(problem; verbose = verbosity > 0)
#     end
#     _, info["model ILP solution time [s]"], _, _, _ = @timed begin
#         optimise!(m, verbosity, timelimit)
#     end
#     _, info["solution extraction time [s]"], _, _, _ = @timed begin
#         design = extract_solution(problem, ind, n1_selected)
#     end
#     info["total time [s]"] = sum([val for (key, val) in info])
#     info["number of variables"] = size(problem)
#     pnull   = problem.toer.score.pnull
#     α       = problem.toer.α
#     mlecomp = mlecompatible(design, pnull, α)
#     if !mlecomp["compatible"] & (verbosity > 0)
#         @warn @sprintf "design is not compatible with MLE-ordering, incompatibility degree is %i/%i" mlecomp["incompatibility degree"] size(mlecomp["details"], 1)
#     end
#     info["MLE-compatible"] = mlecomp
#     return OptimalDesign(design, problem, info)
# end
