abstract type Ordering end

Base.iterate(ordering::Ordering, state = 0) = state > 0 ? nothing : (ordering, state + 1)
Base.length(ordering::Ordering) = 1

Base.show(io::IO, ordering::Ordering) = print(io, string(ordering))
Base.show(io::IO, ::MIME"application/prs.juno.inline", ordering::Ordering) = print(io, string(ordering))



function p_value(x1::TI, x2::TI, p0::TR, ordering::TO, design::TD; orientation::Symbol = :superiority) where {TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign}

    XX = sample_space(design)
    if orientation == :superiority
        inds = larger_or_equal.(XX[:,1], XX[:,2], x1, x2, ordering, design)
    end
    if orientation == :inferiority
        inds = smaller_or_equal.(XX[:,1], XX[:,2], x1, x2, ordering, design)
    end
    return sum(pdf.(XX[:,1], XX[:,2], design, p0)[inds])
end



function compatible(estimator::Estimator, design::AbstractDesign, p0::Real, α::Real)

    !(0 <= p0 <= 1) ? error("p0 must be in [0,1]") : nothing
    XX                      = sample_space(design)
    decisions               = reject_null.(XX[:, 1], XX[:, 2], design)
    inds_reject             = findall(decisions .== 1)
    inds_not_reject         = findall(decisions .== 0)
    ordering                = EstimatorOrdering(estimator)
    pvals                   = p_value.(XX[:,1], XX[:,2], p0, ordering, design; orientation = :superiority)
    max_pval_reject         = pvals[inds_reject] |> maximum
    incompatible_reject     = XX[findall(pvals[inds_reject] .> α), :]
    min_pval_not_reject     = pvals[inds_not_reject] |> minimum
    incompatible_not_rejcet = XX[findall(pvals[inds_not_reject] .<= α), :]
    compatible              = (max_pval_reject < α) & (min_pval_not_reject >= α)
    return compatible, max_pval_reject, min_pval_not_reject, incompatible_reject, incompatible_not_rejcet
end
