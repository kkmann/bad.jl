struct PosteriorMeanEstimator <: Estimator
    prior::Prior
end

string(estimator::PosteriorMeanEstimator) = @sprintf "PosteriorMeanEstimator<%s>" string(estimator.prior)

(estimator::PosteriorMeanEstimator)(x1::Int, x2::Int, design::T) where {T<:AbstractDesign} =
    mean(update(prior, x1 + x2, n(design, x1)))
