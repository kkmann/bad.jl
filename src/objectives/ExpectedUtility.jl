mutable struct ExpectedUtility <: Objective
    prior::Prior
    λ_tp::Real
    λ_fp::Real
    mcr::Real
end

function update!(objective::ExpectedUtility, prior::Prior)
    objective.prior = prior
end

function add!(m, ind, obj::ExpectedUtility, problem::Problem)
    prior_leq_mcr = condition(obj.prior, high = obj.mcr)
    prob_leq_mcr  = cdf(obj.mcr, obj.prior)
    prior_geq_mcr = condition(obj.prior, low = obj.mcr)
    prob_geq_mcr  = 1 - prob_leq_mcr
    ee = @expression(m, sum(
        power(x1, n1, n2, c2, prior_leq_mcr) * dbinom(x1, n1, prior_leq_mcr) * ind[n1, x1, n2, c2] for
            n1 in n1vals(problem),
            x1 in x1vals(n1, problem),
            n2 in n2vals(n1, x1, problem),
            c2 in c2vals(n1, x1, n2, problem)
        )
    )
    ep = @expression(m, sum(
        power(x1, n1, n2, c2, prior_geq_mcr) * dbinom(x1, n1, prior_geq_mcr) * ind[n1, x1, n2, c2] for
            n1 in n1vals(problem),
            x1 in x1vals(n1, problem),
            n2 in n2vals(n1, x1, problem),
            c2 in c2vals(n1, x1, n2, problem)
        )
    )
    ess = @expression(m, sum(
        (n1 + n2) * dbinom(x1, n1, obj.prior) * ind[n1, x1, n2, c2] for
            n1 in n1vals(problem),
            x1 in x1vals(n1, problem),
            n2 in n2vals(n1, x1, problem),
            c2 in c2vals(n1, x1, n2, problem)
        )
    )
    @objective(m, Max,
        obj.λ_fp*prob_leq_mcr*ee + obj.λ_tp*prob_geq_mcr*ep - ess
    )
end
