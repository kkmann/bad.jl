abstract type PowerConstraint end

struct ExpectedPowerConstraint <: PowerConstraint
    threshold::Real
    conditional_threshold::Real
    cprior::Prior
end

minimal_expected_power(prior::Prior, mrv::Real, threshold::Real; conditional_threshold::Real = .5) =
    ExpectedPowerConstraint(threshold, conditional_threshold, condition(prior, low = mrv))

function valid(n1, x1, n2, c2, cnstr::ExpectedPowerConstraint)
    if cnstr.conditional_threshold > 0 # save time if not!
        return !early_stop(c2) & (power(x1, n1, n2, c2, cnstr.cprior) <= cnstr.conditional_threshold) ? false : true
    else
        return true
    end
end


function +(model::DesignIPModel, cnstr::ExpectedPowerConstraint)
    x1vals = model.params["x1vals"]
    n1vals = model.params["n1vals"]
    n2vals = model.params["n2vals"]
    c2vals = model.params["c2vals"]
    valid  = model.params["valid"]
    ind    = model.vars["ind"]
    @constraint(model.jump_model,
        sum(
            power(x1, n1, n2, c2, cnstr.cprior) * predictive_pmf(x1, n1, cnstr.cprior) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        ) >= cnstr.threshold
    )
    return model
end
