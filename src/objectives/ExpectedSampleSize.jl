struct ExpectedSampleSize end
minimize_expected_sample_size() = ExpectedSampleSize()

function +(model::DesignIPModel, objective::ExpectedSampleSize)
    prior = model.params["prior"]
    x1vals = model.params["x1vals"]
    n1vals = model.params["n1vals"]
    n2vals = model.params["n2vals"]
    c2vals = model.params["c2vals"]
    valid  = model.params["valid"]
    ind    = model.vars["ind"]
    @objective(model.jump_model, Min,
        sum(
            (n1 + n2) * predictive_pmf(x1, n1, prior) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        )
    )
    return model
end
