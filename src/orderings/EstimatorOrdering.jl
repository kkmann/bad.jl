struct EstimatorOrdering <: Ordering
    estimator::Estimator
end

string(ordering::EstimatorOrdering) = "EstimatorOrdering"

function strictly_smaller(x1::Int, x2::Int, y1::Int, y2::Int, ordering::EstimatorOrdering, design::T)  where {T<:AbstractDesign}
    ordering.estimator(x1, x2, design) < ordering.estimator(y1, y2, design)
end
function smaller_or_equal(x1::Int, x2::Int, y1::Int, y2::Int, ordering::EstimatorOrdering, design::T)  where {T<:AbstractDesign}
    ordering.estimator(x1, x2, design) <= ordering.estimator(y1, y2, design)
end

function strictly_larger(x1::Int, x2::Int, y1::Int, y2::Int, ordering::EstimatorOrdering, design::T)  where {T<:AbstractDesign}
    ordering.estimator(x1, x2, design) > ordering.estimator(y1, y2, design)
end
function larger_or_equal(x1::Int, x2::Int, y1::Int, y2::Int, ordering::EstimatorOrdering, design::T)  where {T<:AbstractDesign}
    ordering.estimator(x1, x2, design) >= ordering.estimator(y1, y2, design)
end
