struct EstimatorOrdering <: Ordering
    estimator::Estimator
end

string(ordering::EstimatorOrdering) = "EstimatorOrdering"

function (ordering::EstimatorOrdering)(x1::Int, x2::Int, y1::Int, y2::Int, design::T) where {T<:AbstractDesign}
    ordering.estimator(x1, x2, design) < ordering.estimator(y1, y2, design)
end
