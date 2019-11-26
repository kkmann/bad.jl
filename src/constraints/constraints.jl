abstract type TypeOneErrorRateConstraint end

abstract type PowerConstraint end

valid(n1, x1, n2, cnstr::Union{TypeOneErrorRateConstraint, PowerConstraint}) = true
valid(n1, x1, n2, c2, cnstr::Union{TypeOneErrorRateConstraint, PowerConstraint}) = true
