struct MaximumLikelihoodEstimator <: Estimator
end

string(estimator::MaximumLikelihoodEstimator) = "MLE"

function (estimator::MaximumLikelihoodEstimator)(x1::TI, x2::TI, design::TD) where {TI<:Integer,TD<:AbstractDesign}
     return convert(Float64, (x1 + x2) / n(design, x1))
end
