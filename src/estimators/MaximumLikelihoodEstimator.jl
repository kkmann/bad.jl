struct MaximumLikelihoodEstimator <: Estimator
end

string(estimator::MaximumLikelihoodEstimator) = "MaximumLikelihoodEstimator"

function (estimator:: MaximumLikelihoodEstimator)(x1::Int, x2::Union{Int,Nothing}, design::T) where {T<:AbstractDesign}
     return (x1 + x2) / n(design, x1)
end
