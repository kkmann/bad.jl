mutable struct MiniMaxSampleSize <: Objective
    λ::Real
    prior::Prior
end

function (objective::MiniMaxSampleSize)(design::AbstractDesign)

    ess = SampleSize(objective.prior)
    return (1 - objective.λ)*maximum(n(design)) + λ*ess(design)
end

function update(obj::MiniMaxSampleSize, prior::Prior)
    obj = deepcopy(obj)
    obj.prior = prior
    return obj
end

function add!(
        JuMP_model_and_indicator_variables::Tuple,
        objective::MiniMaxSampleSize,
        problem::Problem;
        partial::Tuple{TI,TI} = (0, 0)
    ) where {TI<:Integer}

    m, ind = JuMP_model_and_indicator_variables
    prior  = objective.prior

    @variable(m, nmax)
    for n1_ in problem.n1values
        for x1_ in partial[1]:n1_
            @constraint(m, nmax >= sum(
                (n1_ + n2) * ind[(n1_, x1_, n2, c2)]
                    for (x1, n2, c2) in grid(problem, n1_)
                    if  (x1 == x1_)
                )
            )
        end
    end
    @objective(m, Min,
        (1 - objective.λ)*nmax +
        objective.λ*sum(
            (n1 + n2) * pmf(x1, n1, objective.prior; partial = partial) * ind[(n1, x1, n2, c2)]
                for (n1, x1, n2, c2) in grid(problem)
        )
    )
end
