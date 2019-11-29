mutable struct Problem
    objective::Objective
    toer::TypeOneErrorRateConstraint
    power::PowerConstraint
    n1values:: Vector{Int}
    nmax::Int
    max_rel_increase::Real
    min_rel_increase::Real
    min_abs_increase::Int
    type::Symbol
end

function Problem(
        objective::Objective,
        toer::TypeOneErrorRateConstraint,
        power::PowerConstraint;
        type::Symbol                 = :TwoStage,
        multiple::Real               = 2,
        nmax::Int                    = guess_nmax(p0(toer), α(toer), p1(power), β(power); multiple = multiple),
        n1min_fctr::Real             = .25/multiple,
        n1min::Int                   = type == :OneStage ?  Int(round(nmax/multiple/2)) : Int(round(max(n1min_fctr * nmax, 5))),
        n1max_fctr::Real             = .66/multiple,
        n1max::Int                   = type == :OneStage ? Int(round(nmax)) : Int(round(max(n1min, n1max_fctr * nmax))),
        n1values::Vector{Int}        = zeros(Int, 0),
        max_rel_increase::Real       = 4.,
        min_rel_increase::Real       = 1.1,
        min_abs_increase::Int        = 5
    )
    !any(type .== (:TwoStage, :OneStage, :GroupSequential)) ? error("invalid type specification") : nothing
    nmax > 150 ? (@warn "nmax > 150, consider a continuous approximation using the R package 'adoptr'.") : nothing
    n1values = (length(n1values) > 0) ? n1values : collect(n1min:n1max)
    return Problem(
        objective,
        toer,
        power,
        n1values,
        nmax,
        max_rel_increase,
        min_rel_increase,
        min_abs_increase,
        type
    )
end

function update!(problem::Problem, prior::Prior)
    update!(problem.objective, prior)
    update!(problem.toer, prior)
    update!(problem.power, prior)
end

n1vals(problem::Problem) = problem.n1values
x1vals(n1, problem::Problem) = collect(0:n1)

function n2vals(problem::Problem)
    problem.type == :OneStage ? (return [0]) : nothing
    n1min = min(n1vals(problem))
    n1max = max(n1vals(problem))
    n2min = max(
        problem.min_abs_increase,
        Int(ceil(n1min * problem.min_rel_increase)) - n1min
    )
    n2max = min(
        problem.nmax - n1min,
        Int(floor(n1max * problem.max_rel_increase)) - n1max
    )
    n2 = collect(n2min:n2max)
    return n2min > 0 ? vcat( [0], n2 ) : n2
end
function n2vals(n1, x1, problem::Problem)
    problem.type == :OneStage ? (return [0]) : nothing
    n2min = max(problem.min_abs_increase, Int(ceil(n1*problem.min_rel_increase)) - n1)
    n2max = min(problem.nmax - n1, Int(floor(n1*(problem.max_rel_increase))) - n1)
    n2    = collect(n2min:n2max)
    n2    = n2min > 0 ? vcat( [0], n2 ) : n2
    return n2[valid.(n1, x1, n2, [problem.toer]) .& valid.(n1, x1, n2, [problem.power])]
end

function c2vals(n1, x1, n2, problem::Problem)
    c2 = collect(0:(n2 - 1))
    c2 = vcat( [-Inf], c2, [Inf])
    return c2[valid.(n1, x1, n2, c2, [problem.toer]) .& valid.(n1, x1, n2, c2, [problem.power])]
end

function get_IP_model(problem::Problem)
    prog = ProgressMeter.Progress(
        1 + length(n1vals(problem)) + 4,
        desc = "Building IP Problem: ",
        dt = 0.5,
        barglyphs = ProgressMeter.BarGlyphs("[=> ]"),
        barlen = 30,
        color = :gray
    )
    ProgressMeter.next!(prog; showvalues = [(Symbol("building model"), "...")], valuecolor = :gray)
    m = Model()
    @variable(m, # ind[x1, n1, n2, c2] == 1 iff n1 = n1, n2(x1) = n2, c2(x1) = c2
        ind[
            n1 in n1vals(problem),
            x1 in x1vals(n1, problem),
            n2 in n2vals(n1, x1, problem),
            c2 in c2vals(n1, x1, n2, problem)
        ],
        Bin
    )
    # need to make sure that exactly one n1 value is selected, IP or trick via
    # auxiliary variables: n1_selected[n1] == 1 iff n1 = n1
    @variable(m, n1_selected[n1 in n1vals(problem)], Bin)
    @constraint(m, sum(n1_selected[n1] for n1 in n1vals(problem)) == 1)
    # implement all constraints conditional on n1
    for n1 in n1vals(problem)
        ProgressMeter.next!(prog; showvalues = [(Symbol("building constraints for n1"), n1)], valuecolor = :gray)
        # make sure that n1_selected[n1] is 1 iff any of the other values is assigned
        @constraint(m,
            2*problem.nmax * n1_selected[n1] >= sum(
                ind[n1, x1, n2, c2] for
                x1 in x1vals(n1, problem), n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
            )
        )
        for x1 in x1vals(n1, problem)
            # make its a function in x1
            @constraint(m,
                n1_selected[n1] == sum(
                    ind[n1, x1, n2, c2] for
                    n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
                )
            )
            # make sure that n1(x1) is constant
            if 0 < x1 <= n1
                @constraint(m,
                    sum(ind[n1,     x1, n2, c2] for
                        n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
                    ) ==
                    sum(ind[n1, x1 - 1, n2, c2] for
                        n2 in n2vals(n1, x1 - 1, problem), c2 in c2vals(n1, x1 - 1, n2, problem)
                    )
                )
                if problem.type == :GroupSequential
                    for n2 in n2vals(n1, x1, problem)
                        @constraint(m,
                            sum( ind[n1, x1, n2, c2] for c2 in c2vals(n1, x1, n2, problem) if isfinite(c2) ) +
                            sum( ind[n1, x1, nn2, c2] for nn2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, nn2, problem) if c2 == -Inf )>=
                            sum( ind[n1, x1 - 1, n2, c2] for c2 in c2vals(n1, x1 - 1, n2, problem) if isfinite(c2) )
                        )
                    end
                end
            end
            # make sure we have contiguous stopping for efficacy / futility
            if x1 > 0
                @constraint(m,
                    ind[n1, x1 - 1, 0, Inf] >= ind[n1, x1, 0, Inf]
                )
            end
            if x1 < n1
                @constraint(m,
                    ind[n1, x1 + 1, 0, -Inf] >= ind[n1, x1, 0, -Inf]
                )
            end
        end
    end
    ProgressMeter.next!(prog; showvalues = [(Symbol("adding constraint"), "type one error rate")], valuecolor = :gray)
    add!(m, ind, problem.toer, problem)
    ProgressMeter.next!(prog; showvalues = [(Symbol("adding constraint"), "power")], valuecolor = :gray)
    add!(m, ind, problem.power, problem)
    ProgressMeter.next!(prog; showvalues = [(Symbol("adding objective"), "...")], valuecolor = :gray)
    add!(m, ind, problem.objective, problem)
    ProgressMeter.next!(prog; showvalues = [(Symbol("finishing up"), "...done!")], valuecolor = :gray)
    return m, ind, n1_selected
end



mutable struct OptimalDesign <: AbstractDesign
    n2::Vector{Int}
    c2::Vector{Real}
    model::Problem
    score::Real
    info::Dict{String,Any}
end



function optimise(problem::Problem; verbosity = 3, timelimit = 300)
    tick()
    m, ind, n1_selected = get_IP_model(problem)
    time_problem_generation = tok()
    tick()
    optimize!(
        m,
        with_optimizer(GLPK.Optimizer, msg_lev = verbosity, tm_lim = 1000*timelimit)
    )
    time_problem_solution = tok()
    # termination_status(model.jump_model)
    vals = value.(ind)
    # find n1
    n1 = n1vals(problem)[findfirst(value.(n1_selected).data .== 1.0)]
    # extract solution
    c2_res = repeat([Inf], n1 + 1)
    n2_res = zeros(n1 + 1)
    for x1 in x1vals(n1, problem)
        cntr = 0 # make sure we have a proper solutipon
        for n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
            if vals[n1, x1, n2, c2] == 1
                c2_res[x1 + 1] = c2
                n2_res[x1 + 1] = n2
                cntr          += 1
            end
        end
        cntr != 1 ? error() : nothing
    end
    score = objective_value(m)
    info  = Dict(
        "#variables" => length(ind),
        "problem generation time" => time_problem_generation,
        "solution time" => time_problem_solution
    )
    return OptimalDesign(n2_res, c2_res, problem, score, info)
end
