mutable struct MiniMaxSampleSize <: Objective
    λ::Real
    prior::Prior
end

function (objective::MiniMaxSampleSize)(design::AbstractDesign)

    return maximum(n(design)) + λ*expectation(p -> expected_sample_size(desing, p), objective.prior)
end

update!(objective::MiniMaxSampleSize, prior::Prior) = nothing

function add!(jump_model, ind, objective::MiniMaxSampleSize, problem::Problem)
    @variable(jump_model, nmax >= 0)
    for n1 in n1vals(problem)
        for x1 in 0:n1
            @constraint(jump_model, nmax >= sum(
                (n1 + n2) * ind[n1, x1, n2, c2] for
                    n2 in n2vals(n1, x1, problem),
                    c2 in c2vals(n1, x1, n2, problem)
                )
            )
        end
    end
    @objective(jump_model, Min,
        (1 - objective.λ)*nmax + objective.λ*sum(
            (n1 + n2) * dbinom(x1, n1, objective.prior) * ind[n1, x1, n2, c2] for
                n1 in n1vals(problem),
                x1 in x1vals(n1, problem),
                n2 in n2vals(n1, x1, problem),
                c2 in c2vals(n1, x1, n2, problem)
        )
    )
end
