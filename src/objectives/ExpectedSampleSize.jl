struct ExpectedSampleSize <: Objective
    prior::Prior
end
minimise_expected_sample_size(prior) = ExpectedSampleSize(prior)

function add!(jump_model, ind, objective::ExpectedSampleSize, problem::Problem)
    @objective(jump_model, Min,
        sum(
            (n1 + n2) * dbinom(x1, n1, objective.prior) * ind[n1, x1, n2, c2] for
                n1 in n1vals(problem),
                x1 in x1vals(n1, problem),
                n2 in n2vals(n1, x1, problem),
                c2 in c2vals(n1, x1, n2, problem)
        )
    )
end
