struct PosteriorMean <: Estimator
    prior::Prior
end

string(estimator::PosteriorMean) = @sprintf "PosteriorMean<%s>" string(estimator.prior)

function (estimator::PosteriorMean)(x1::TI, x2::TI, design::TD) where {TI<:Integer,TD<:AbstractDesign}
    mean(update(estimator.prior, x1 + x2, n(design, x1)))
end


# TODO: implement this as cached() / CachedEstimator for all!
struct PosteriorMeanPrecalculated{TI<:Integer,TR<:Real,TD<:AbstractDesign,TP<:Prior} <: Estimator
    design::TD
    prior::TP
    XX::Array{TI,2} # sample space
    estimates::Array{TR,1}
end

string(estimator::PosteriorMeanPrecalculated) = @sprintf "PosteriorMean(pre-calculated)<%s;%s>" string(estimator.prior) string(estimator.design)

function (estimator::PosteriorMeanPrecalculated)(x1::TI, x2::TI, design::TD) where {TI<:Integer,TD<:AbstractDesign}

    for i in 1:size(estimator.XX, 1)
        if [x1, x2] == estimator.XX[i,:]
            return estimator.estimates[i]
        end
    end
    error("(x1,x2) not found in sample space, valid observation?")
end

function PosteriorMeanPrecalculated(prior::Prior, design::AbstractDesign)

    XX        = sample_space(design)
    estimates = zeros(size(XX, 1))
    for i in 1:size(XX, 1)
        x1, x2 = XX[i,:]
        estimates[i] = mean(update(prior, x1 + x2, n(design, x1)))
    end
    return PosteriorMeanPrecalculated{eltype(XX),eltype(estimates),typeof(design),typeof(prior)}(
        design, prior, XX, estimates
    )
end
