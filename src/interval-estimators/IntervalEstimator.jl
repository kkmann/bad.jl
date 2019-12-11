abstract type IntervalEstimator end

Base.iterate(estimator::IntervalEstimator, state = 0) = state > 0 ? nothing : (estimator, state + 1)
Base.length(estimator::IntervalEstimator) = 1

Base.show(io::IO, estimator::IntervalEstimator) = print(io, string(estimator))
Base.show(io::IO, ::MIME"application/prs.juno.inline", estimator::IntervalEstimator) = print(io, string(estimator))

Base.string(estimator::IntervalEstimator) = string(typeof(estimator))



function get_bounds(estimator::IntervalEstimator, x1::TI, x2::TI) where {TI<:Integer}

    for i in 1:size(estimator.XX, 1)
        if [x1, x2] == estimator.XX[i,:]
            return estimator.bounds[i, :]
        end
    end
    error("(x1,x2) not found in sample space, valid observation?")
end
function get_bounds(estimator::IntervalEstimator, x1::Vector{TI}, x2::Vector{TI}) where {TI<:Integer}

    return convert(Matrix{eltype(estimator.bounds)}, hcat(get_bounds.(estimator, x1, x2)...)')
end

(estimator::IntervalEstimator)(x1::TI, x2::TI) where {TI<:Integer} = get_bounds(estimator, x1, x2)
(estimator::IntervalEstimator)(x1::Vector{TI}, x2::Vector{TI}) where {TI<:Integer} = get_bounds(estimator, x1, x2)



function compatible(estimator::IntervalEstimator, design::AbstractDesign, pnull::Real)

    @assert (0 <= pnull <= 1) "pnull must be in [0,1]"

    XX                      = sample_space(design)
    design_rejects          = reject.(XX[:,1], XX[:,2], design)
    lower_boundaries        = get_bounds(estimator, XX[:,1], XX[:,2])[:,1]
    iv_rejects              = lower_boundaries .> pnull
    incompatibility_degree  = sum(design_rejects .& .!iv_rejects)

    df = DataFrames.DataFrame(
        x1 = XX[:,1],
        x2 = XX[:,2],
        design_rejects = design_rejects,
        ci_reject = iv_rejects,
        lower_boundaries = lower_boundaries
    )

    return Dict{String,Any}(
        "compatible" => incompatibility_degree == 0,
        "incompatibility degree" => incompatibility_degree,
        "details" => df
    )
end


function coverage_probability(estimator::IntervalEstimator, p::Real)

    XX    = estimator.XX
    pmff  = pmf.(XX[:,1], XX[:,2], estimator.design, p)
    lower = sum(pmff[estimator.bounds[:,1] .<= p])
    upper = sum(pmff[estimator.bounds[:,2] .>= p])
    joint = sum(pmff[(estimator.bounds[:,1] .<= p) .& (estimator.bounds[:,2] .>= p)])
    return [lower, joint, upper]
end

function coverage_probability(estimator::IntervalEstimator, p::Vector{TR}) where {TR<:Real}

    return convert(Matrix{TR}, hcat(coverage_probability.(estimator, p)...)')
end


function mean_width(estimator::IntervalEstimator, p::Real)

    pmff   = pmf.(estimator.XX[:,1], estimator.XX[:,2], estimator.design, p)
    return sum( (estimator.bounds[:,2] .- estimator.bounds[:,1]) .* pmff )
end
