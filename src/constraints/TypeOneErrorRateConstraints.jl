abstract type TypeOneErrorRateConstraint end

mutable struct MaximalTypeOneErrorRateConstraint{T<:Real}
    p0::T
    threshold::T
end

maximal_type_one_error_rate(p0, α) = MaximalTypeOneErrorRateConstraint(p0, α)

function valid(n1, x1, n2, c2, cnstr::MaximalTypeOneErrorRateConstraint)
    if n2 == 0
        return x1 > findfirst(1 .- pbinom.(0:n1, n1, p0) .<= α) + k ? false : true
    else
        return true
    end
end

function +(model::DesignIPModel, cnstr::MaximalTypeOneErrorRateConstraint)
    x1vals = model.params["x1vals"]
    n1vals = model.params["n1vals"]
    n2vals = model.params["n2vals"]
    c2vals = model.params["c2vals"]
    valid  = model.params["valid"]
    ind    = model.vars["ind"]
    @constraint(model.jump_model,
        sum(
            power(x1, n1, n2, c2, cnstr.p0) * predictive_pmf(x1, n1, cnstr.p0) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        ) <= cnstr.threshold
    )
    return model
end
