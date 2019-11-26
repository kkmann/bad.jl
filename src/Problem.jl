struct Problem
    objective::Objective
    toer::TypeOneErrorRateConstraint
    power::PowerConstraint
    params::Dict{String,Any}
end

function Problem(
        objective::Objective,
        toer::TypeOneErrorRateConstraint,
        power::PowerConstraint;
        type::Symbol                 = :TwoStage,
        multiple::Real               = type == :OneStage ? 1.25 : 2,
        nmax::Int                    = guess_nmax(p0(toer), α(toer), p1(power), β(power); multiple = multiple),
        n1min_fctr::Real             = .25/multiple,
        n1min::Int                   = type == :OneStage ? Int(round(max(.66*nmax, 5))) : Int(round(max(n1min_fctr * nmax, 5))),
        n1max_fctr::Real             = .66/multiple,
        n1max::Int                   = type == :OneStage ? nmax : Int(round(max(n1min, n1max_fctr * nmax))),
        n1values::Vector{Int}        = zeros(Int, 0),
        max_rel_increase::Real       = 3.,
        min_rel_increase::Real       = 1.1,
        min_abs_increase::Int        = 5
    )

    !any(type .== (:TwoStage, :OneStage, :GroupSequential)) ? error("invalid type specification") : nothing
    nmax > 200 ? (@warn "nmax > 200, consider a continuous approximation using the R package 'adoptr'.") : nothing

    return Problem(
        objective,
        toer,
        power,
        Dict(
            "type"             => type,
            "nmax"             => nmax,
            "n1min"            => n1min,
            "n1max"            => n1max,
            "n1values"         => n1values,
            "max_rel_increase" => max_rel_increase,
            "min_rel_increase" => min_rel_increase,
            "min_abs_increase" => min_abs_increase
        )
    )
end

n1vals(problem::Problem) = (length(problem.params["n1values"]) > 0) ? problem.params["n1values"] : collect(problem.params["n1min"]:problem.params["n1max"])
x1vals(n1, problem::Problem) = collect(0:n1)

function n2vals(problem::Problem)
    n2min = max(
        problem.params["min_abs_increase"],
        Int(ceil(n1*problem.params["min_rel_increase"])) - problem.params["n1min"]
    )
    n2max = min(
        problem.params["nmax"] - problem.params["n1min"],
        Int(floor(n1 * (problem.params["max_rel_increase"]))) - problem.params["n1min"]
    )
    return n2min > 0 ? cat( [0], collect(n2min:n2max) ) : collect(n2min:n2max)
end
function n2vals(n1, x1, problem::Problem)
    n2min = max(problem.params["min_abs_increase"], Int(ceil(n1*problem.params["min_rel_increase"])) - n1)
    n2max = min(problem.params["nmax"] - n1, Int(floor(n1 * (problem.params["max_rel_increase"]))) - n1)
    n2 = collect(n2min:n2max)
    n2    = n2min > 0 ? vcat( [0], n2 ) : n2
    return n2[valid.(n1, x1, n2, [problem.toer]) .& valid.(n1, x1, n2, [problem.power])]
end

function c2vals(n1, x1, n2, problem::Problem)
    c2 = collect(0:(n2 - 1))
    c2 = vcat( [EarlyFutility, EarlyEfficacy], c2 )
    return c2[valid.(n1, x1, n2, c2, [problem.toer]) .& valid.(n1, x1, n2, c2, [problem.power])]
end

function get_IP_model(problem::Problem)
    # n1vals, x1vals, n2vals, c2vals, c2vals_stop, c2vals_cont = get_grid(problem)
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
    if problem.params["type"] == :GroupSequential # same for n2 in group sequential case
        @variable(m, n2_selected[n2 in n2vals(problem)], Bin)
        @constraint(m, sum(n2_selected[n2] for n2 in n2vals(problem)) == 1)
    end
    # implement all constraints conditional on n1
    @showprogress for n1 in n1vals(problem)
        # make sure that n1_selected[n1] is 1 iff any of the other values is assigned
        @constraint(m,
            2*problem.params["nmax"] * n1_selected[n1] >= sum(
                ind[n1, x1, n2, c2] for
                x1 in x1vals(n1, problem), n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
            )
        )
        # make sure that n2_selected[n2] is 1 iff any of the other values is assigned
        if problem.params["type"] == :GroupSequential
            for n2 in n2vals
                @constraint(m,
                    2*problem.params["nmax"] * n2_selected[n2] >= sum(
                        ind[n1, x1, n2, c2] for
                        x1 in x1vals(n1, problem), n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
                    )
                )
            end
        end
        for x1 in x1vals(n1, problem)
            # make its a function in x1
            @constraint(m,
                n1_selected[n1] == sum(
                    ind[n1, x1, n2, c2] for
                    n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
                )
            )
            # make sure that n1(x1) is constant
            if 0 < x1 < n1
                @constraint(m,
                    sum(ind[n1,     x1, n2, c2] for
                        n2 in n2vals(n1, x1, problem), c2 in c2vals(n1, x1, n2, problem)
                    ) ==
                    sum(ind[n1, x1 - 1, n2, c2] for
                        n2 in n2vals(n1, x1 - 1, problem), c2 in c2vals(n1, x1 - 1, n2, problem)
                    )
                )
            end
            # make sure we have contiguous stopping for efficacy / futility
            if x1 > 0
                @constraint(m,
                    ind[n1, x1 - 1, 0, EarlyFutility] >= ind[n1, x1, 0, EarlyFutility]
                )
            end
            if x1 < n1
                @constraint(m,
                    ind[n1, x1 + 1, 0, EarlyEfficacy] >= ind[n1, x1, 0, EarlyEfficacy]
                )
            end
        end
    end
    add!(m, ind, problem.toer, problem)
    add!(m, ind, problem.power, problem)
    add!(m, ind, problem.objective, problem)
    return m, ind, n1_selected
end



mutable struct OptimalDesign <: AbstractDesign
    n2::Vector{Int}
    c2::Vector{CriticalValue}
    model::Problem
    score::Real
    # check n2/c2
end



function optimise(problem::Problem; verbosity = 3, timelimit = 300)
    # n1vals, x1vals, n2vals, c2vals, c2vals_stop, c2vals_cont = get_grid(problem)
    m, ind, n1_selected = get_IP_model(problem)
    optimize!(
        m,
        with_optimizer(GLPK.Optimizer, msg_lev = verbosity, tm_lim = 1000*timelimit)
    )
    # termination_status(model.jump_model)
    vals = value.(ind)
    # find n1
    n1 = n1vals(problem)[findfirst(value.(n1_selected).data .== 1.0)]
    # extract solution
    c2_res = convert(Vector{CriticalValue}, repeat([EarlyFutility], n1 + 1))
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
    return OptimalDesign(n2_res, c2_res, problem, score)
end
