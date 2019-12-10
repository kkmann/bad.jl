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

function more_extreme(x1::TI, x2::TI, x1_::TI, x2_::TI,
    ordering::EstimatorOrdering, design::TD; orientation = nothing)  where {TI<:Integer,TD<:AbstractDesign}

    orientation = (orientation == nothing) ? ordering.orientation : orientation
    if orientation == :superiority
        return ordering.estimator(x1, x2, design) >= ordering.estimator(x1_, x2_, design)
    end
    if orientation == :inferiority
        return ordering.estimator(x1, x2, design) <= ordering.estimator(x1_, x2_, design)
    end
    error("orientation must be :superiororit or :inferiority")
end

function more_extreme(x1::TI, x2::TI, x1_::TI, x2_::TI, estimator::TE,
    design::TD; orientation = nothing)  where {TI<:Integer,TE<:Estimator,TD<:AbstractDesign}

    return more_extreme(x1, x2, x1_, x2_, EstimatorOrdering(estimator), design,
        orientation = oreintation)
end



function p_value(x1::TI, x2::TI, p0::TR, ordering::TO, design::TD; orientation = nothing) where {TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign}

    XX   = sample_space(design)
    inds = more_extreme.(XX[:,1], XX[:,2], x1, x2, ordering, design, orientation = orientation)
    return min(1, max(0, sum(pdf.(XX[:,1], XX[:,2], design, p0)[inds]) ) )
end

function p_value(x1::TI, x2::TI, p0::TR, estimator::TE, design::TD; orientation::Symbol = :superiority) where {TI<:Integer,TR<:Real,TE<:Estimator,TD<:AbstractDesign}

    return p_value(x1, x2, p0, EstimatorOrdering(estimator, orientation), design)
end



function compatible(ordering::Ordering, design::AbstractDesign, pnull::Real, α::Real)

    @assert (0 <= pnull <= 1) "pnull must be in [0,1]"

    XX                      = sample_space(design)
    design_rejects          = reject.(XX[:,1], XX[:,2], design)
    pvals                   = p_value.(XX[:,1], XX[:,2], pnull, ordering, design)
    ordering_rejects        = pvals .< α
    incompatibility_degree  = sum(design_rejects .& .!ordering_rejects)

    df = DataFrames.DataFrame(
        x1 = XX[:,1],
        x2 = XX[:,2],
        design_rejects = design_rejects,
        ordering_rejects = ordering_rejects,
        p_value = pvals
    )

    return Dict{String,Any}(
        "compatible" => incompatibility_degree == 0,
        "incompatibility degree" => incompatibility_degree,
        "details" => df
    )
end

function mlecompatible(design::AbstractDesign, pnull::Real, α::Real)

    return compatible(EstimatorOrdering(MaximumLikelihoodEstimator()), design, pnull, α)
end
