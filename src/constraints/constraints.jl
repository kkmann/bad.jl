abstract type Constraint end

abstract type TypeOneErrorRateConstraint <: Constraint end

p0(cnstr::TypeOneErrorRateConstraint) = error("not implemented")
α(cnstr::TypeOneErrorRateConstraint) = error("not implemented")

abstract type PowerConstraint <: Constraint end

p1(cnstr::PowerConstraint) = error("not implemented")
β(cnstr::PowerConstraint) = error("not implemented")

(::Constraint)(design::AbstractDesign) = error("not implemented")
(::Constraint)(design::AbstractDesign, x1::Integer) = error("not implemented")

valid(n1, x1, n2, cnstr::Constraint) = true
valid(n1, x1, n2, c2, cnstr::Constraint) = true

update!(cnstr::Constraint, prior::Prior) = cnstr

# add!(m, ind, cnstr::Constraint, problem::Problem) = nothing
# add!(m, ind, cnstr::Constraint, problem::Problem, xx1, nn1, old_design::AbstractDesign) = nothing
