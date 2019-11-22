struct ExpectedPowerConstraint
end
expected_power_constraint() = ExpectedPowerConstraint()

function +(model::DesignIPModel, cnstr::ExpectedPowerConstraint)
    β      = model.params["β"]
    cprior = model.params["cprior"]
    x1vals = model.params["x1vals"]
    n1vals = model.params["n1vals"]
    n2vals = model.params["n2vals"]
    c2vals = model.params["c2vals"]
    valid  = model.params["valid"]
    ind    = model.vars["ind"]
    @constraint(model.jump_model,
        sum(
            power(x1, n1, n2, c2, cprior) * predictive_pmf(x1, n1, cprior) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        ) >= 1 - β
    )
    return model
end
