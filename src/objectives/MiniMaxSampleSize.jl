mutable struct MiniMaxSampleSize <: Objective
end

(objective::MiniMaxSampleSize)(design::AbstractDesign) = maximum(n(design))

update!(objective::MiniMaxSampleSize, prior::Prior) = nothing

function add!(jump_model, ind, objective::MiniMaxSampleSize, problem::Problem)
    @variable(jump_model, nmax >= 0)
    @constraint(jump_model, nmax >= sum(
        (n1 + n2) * ind[n1, x1, n2, c2] for
            n1 in n1vals(problem),
            x1 in x1vals(n1, problem),
            n2 in n2vals(n1, x1, problem),
            c2 in c2vals(n1, x1, n2, problem)
        )
    )
    @objective(jump_model, Min, nmax)
end
