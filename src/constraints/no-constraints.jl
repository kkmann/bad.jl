add!(jump_model, ind, cnstr::Union{PowerConstraint,TypeOneErrorRateConstraint}, problem::Problem) = nothing
add!(jump_model, ind, cnstr::Union{PowerConstraint,TypeOneErrorRateConstraint}, problem::Problem, xx1, nn1, old_design::AbstractDesign) =
    nothing

struct NoPowerConstraint <: PowerConstraint
    p1::Real
    β::Real
end
p1(cnstr::NoPowerConstraint) = cnstr.p1
β(cnstr::NoPowerConstraint) = cnstr.β
function valid(n1, x1, n2, cnstr::NoPowerConstraint)
    ( (power(x1, n1, cnstr.p1) > 1 - cnstr.β) & (n2 > 0) ) ?
            (return false) : true
end
function valid(n1, x1, n2, c2, cnstr::NoPowerConstraint)
    if n2 > 0 # not stopping early
        conditional_power = power(n2, c2, cnstr.p1)
        conditional_power <= (1 - cnstr.β) / 2 ? (return false) : nothing
        conditional_power >= .99 ?  (return false) : nothing
    end
    return true
end

struct NoTypeOneErrorRateConstraint <: TypeOneErrorRateConstraint
    p0::Real
    α::Real
end
p0(cnstr::NoTypeOneErrorRateConstraint) = cnstr.p0
α(cnstr::NoTypeOneErrorRateConstraint) = cnstr.α
function valid(n1, x1, n2, cnstr::NoTypeOneErrorRateConstraint)
    # if we continue, first stage p value must be rather large (buffer for multiple testing)
    ( (x1 > findfirst(1 .- pbinom.(0:n1, n1, cnstr.p0) .<= cnstr.α) + 3) & (n2 > 0) ) ?
        (return false) : nothing
    return true
end
