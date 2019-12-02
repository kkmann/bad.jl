struct PosteriorMean <: Estimator
    prior::Prior
end

string(estimator::PosteriorMean) = @sprintf "PosteriorMean<%s>" string(estimator.prior)

function (estimator::PosteriorMean)(x1::TI, x2::TI, design::TD) where {TI<:Integer,TD<:AbstractDesign}
    mean(update(estimator.prior, x1 + x2, n(design, x1)))
end
