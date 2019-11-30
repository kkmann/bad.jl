abstract type Constraint end

abstract type TypeOneErrorRateConstraint <: Constraint end

abstract type PowerConstraint <: Constraint end

valid(n1, x1, n2, cnstr::Union{TypeOneErrorRateConstraint, PowerConstraint}) = true
valid(n1, x1, n2, c2, cnstr::Union{TypeOneErrorRateConstraint, PowerConstraint}) = true

update!(cnstr::Union{TypeOneErrorRateConstraint, PowerConstraint}, prior::Prior) = cnstr
