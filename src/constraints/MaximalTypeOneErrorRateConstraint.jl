mutable struct MaximalTypeOneErrorRateConstraint <: TypeOneErrorRateConstraint
    p0::Real
    threshold::Real
    k::Int
end

maximal_type_one_error_rate(p0, α; k::Int = 2) =
    MaximalTypeOneErrorRateConstraint(p0, α, k)

function valid(n1, x1, n2, cnstr::MaximalTypeOneErrorRateConstraint)
    # if we continue, first stage p value must be rather large (buffer for multiple testing)
    ( (x1 > findfirst(1 .- pbinom.(0:n1, n1, cnstr.p0) .<= cnstr.threshold) + cnstr.k) & (n2 > 0) ) ?
        (return false) : nothing
    return true
end

p0(cnstr::MaximalTypeOneErrorRateConstraint) = cnstr.p0
α(cnstr::MaximalTypeOneErrorRateConstraint) = cnstr.threshold

function add!(jump_model, ind, cnstr::MaximalTypeOneErrorRateConstraint, problem::Problem)
    @constraint(jump_model,
        sum(
            power(n2, c2, cnstr.p0) * dbinom(x1, n1, cnstr.p0) * ind[n1, x1, n2, c2] for
                n1 in n1vals(problem),
                x1 in x1vals(n1, problem),
                n2 in n2vals(n1, x1, problem),
                c2 in c2vals(n1, x1, n2, problem)
        ) <= cnstr.threshold
    )
end
