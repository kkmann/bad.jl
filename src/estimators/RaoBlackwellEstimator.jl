struct RaoBlackwellEstimator <: Estimator
end

string(estimator::RaoBlackwellEstimator) = "RaoBlackwellEstimator"

function (estimator::RaoBlackwellEstimator)(x1::TI, x2::TI, design::TD) where {TI<:Integer,TD<:AbstractDesign}
    x = x1 + x2
    xx1   = 0:n1(design)
    xx1   = xx1[n2.(design, xx1) .== n2(design, 1)]
    nom   = sum( binomial.(n1(design) .- 1, xx1 .- 1) .* binomial.(n2.(design, xx1), x .- xx1) )
    denom = sum( binomial.(n1(design), xx1) .* binomial.(n2.(design, xx1), x - xx1) )
    return nom / denom
end
