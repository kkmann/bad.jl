mutable struct ExpectedPowerConstraint <: PowerConstraint
    threshold::Real
    conditional_threshold::Real
    cprior::Prior
    power_curtail::Real
end

minimal_expected_power(prior::Prior, mrv::Real, threshold::Real;
    conditional_threshold::Real = .5, power_curtail = .999) =
    ExpectedPowerConstraint(threshold, conditional_threshold, condition(prior, low = mrv), power_curtail)

function (cnstr::ExpectedPowerConstraint)(design::AbstractDesign)
    power(design, cnstr.cprior) - cnstr.threshold
end
function (cnstr::ExpectedPowerConstraint)(design::AbstractDesign, x1::Integer)
    power(x1, design, cnstr.cprior) - cnstr.conditional_threshold
end

function valid(n1, x1, n2, cnstr::ExpectedPowerConstraint)
    ( (power(x1, n1, cnstr.cprior) > cnstr.threshold) & (n2 > 0) ) ?
            (return false) : true
end

function valid(n1, x1, n2, c2, cnstr::ExpectedPowerConstraint)
    if n2 > 0 # not stopping early
        conditional_power = power(x1, n1, n2, c2, cnstr.cprior)
        conditional_power <= cnstr.conditional_threshold ? (return false) : nothing
        conditional_power >= cnstr.power_curtail ?  (return false) : nothing
    end
    return true
end

p1(cnstr::ExpectedPowerConstraint) = mean(cnstr.cprior)
β(cnstr::ExpectedPowerConstraint) = 1 - cnstr.threshold

function update!(cnstr::ExpectedPowerConstraint, prior::Prior; mrv = cnstr.cprior.low)
    cnstr.cprior = condition(prior, low = mrv)
end

function add!(jump_model, ind, cnstr::ExpectedPowerConstraint, problem::Problem)
    @constraint(jump_model,
        sum(
            power(x1, n1, n2, c2, cnstr.cprior) * dbinom(x1, n1, cnstr.cprior) * ind[n1, x1, n2, c2] for
                n1 in n1vals(problem),
                x1 in x1vals(n1, problem),
                n2 in n2vals(n1, x1, problem),
                c2 in c2vals(n1, x1, n2, problem)
        ) >= cnstr.threshold
    )
end

function add!(jump_model, ind, cnstr::ExpectedPowerConstraint, problem::Problem, xx1, nn1, old_design::AbstractDesign)
    ccprior = update(cnstr.cprior, xx1, nn1)
    @constraint(jump_model,
        sum(
            power(x1, n1, n2, c2, ccprior) * dbinom(x1 - xx1, n1 - nn1, ccprior) * ind[n1, x1, n2, c2] for
                n1 in n1vals(problem),
                x1 in x1vals(n1, problem),
                n2 in n2vals(n1, x1, problem),
                c2 in c2vals(n1, x1, n2, problem) if
            n1 >= nn1
        ) >= power(cnstr.cprior, xx1, nn1, old_design)
    )
end
