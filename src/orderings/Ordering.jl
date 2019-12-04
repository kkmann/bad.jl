abstract type Ordering end

Base.iterate(ordering::Ordering, state = 0) = state > 0 ? nothing : (ordering, state + 1)
Base.length(ordering::Ordering) = 1

Base.show(io::IO, ordering::Ordering) = print(io, string(ordering))
Base.show(io::IO, ::MIME"application/prs.juno.inline", ordering::Ordering) = print(io, string(ordering))



struct EstimatorOrdering{TE<:Estimator} <: Ordering
    estimator::TE
    orientation::Symbol
end

function EstimatorOrdering(estimator::Estimator; orientation::Symbol = :superiority)
    EstimatorOrdering{typeof(estimator)}(estimator, orientation)
end

string(ordering::EstimatorOrdering) = "EstimatorOrdering"

function more_extreme(x1::TI, x2::TI, x1_::TI, x2_::TI, ordering::EstimatorOrdering, design::TD)  where {TI<:Integer,TD<:AbstractDesign}

    if ordering.orientation == :superiority
        return ordering.estimator(x1, x2, design) >= ordering.estimator(x1_, x2_, design)
    end
    if ordering.orientation == :inferiority
        return ordering.estimator(x1, x2, design) <= ordering.estimator(x1_, x2_, design)
    end
    error("orientation must be :superiororit or :inferiority")
end

function more_extreme(x1::TI, x2::TI, x1_::TI, x2_::TI, estimator::TE, design::TD)  where {TI<:Integer,TE<:Estimator,TD<:AbstractDesign}

    return more_extreme(x1, x2, x1_, x2_, EstimatorOrdering(estimator), design)
end



function p_value(x1::TI, x2::TI, p0::TR, ordering::TO, design::TD) where {TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign}

    XX   = sample_space(design)
    inds = more_extreme.(XX[:,1], XX[:,2], x1, x2, ordering, design)
    return min(1, max(0, sum(pdf.(XX[:,1], XX[:,2], design, p0)[inds]) ) )
end

function p_value(x1::TI, x2::TI, p0::TR, estimator::TE, design::TD; orientation::Symbol = :superiority) where {TI<:Integer,TR<:Real,TE<:Estimator,TD<:AbstractDesign}

    XX   = sample_space(design)
    inds = more_extreme.(XX[:,1], XX[:,2], x1, x2, EstimatorOrdering(estimator, orientation), design)
    return min(1, max(0, sum(pdf.(XX[:,1], XX[:,2], design, p0)[inds]) ) )
end



function compatible(estimator::Estimator, design::AbstractDesign, p0::Real, α::Real)

    !(0 <= p0 <= 1) ? error("p0 must be in [0,1]") : nothing
    XX                      = sample_space(design)
    decisions               = reject_null.(XX[:, 1], XX[:, 2], design)
    inds_reject             = findall(decisions .== 1)
    inds_not_reject         = findall(decisions .== 0)
    ordering                = EstimatorOrdering(estimator; orientation = :superiority)
    pvals                   = p_value.(XX[:,1], XX[:,2], p0, ordering, design)
    max_pval_reject         = pvals[inds_reject] |> maximum
    incompatible_reject     = XX[findall(pvals[inds_reject] .> α), :]
    min_pval_not_reject     = pvals[inds_not_reject] |> minimum
    incompatible_not_rejcet = XX[findall(pvals[inds_not_reject] .<= α), :]
    compatible              = (max_pval_reject < α) & (min_pval_not_reject >= α)
    return compatible, max_pval_reject, min_pval_not_reject, incompatible_reject, incompatible_not_rejcet
end
