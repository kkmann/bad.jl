abstract type IntervalEstimator end

Base.iterate(estimator::IntervalEstimator, state = 0) = state > 0 ? nothing : (ci, state + 1)
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
(estimator::IntervalEstimator)(x1::TI, x2::TI) where {TI<:Integer} = get_bounds(estimator, x1, x2)



function coverage_probability(estimator::IntervalEstimator, p::Real)

    XX    = estimator.XX
    pmf   = pdf.(XX[:,1], XX[:,2], estimator.design, p)
    lower = sum(pmf[estimator.bounds[:,1] .<= p])
    upper = sum(pmf[estimator.bounds[:,2] .>= p])
    joint = sum(pmf[(estimator.bounds[:,1] .<= p) .& (estimator.bounds[:,2] .>= p)])
    return [lower, joint, upper]
end



function mean_width(estimator::IntervalEstimator, p::Real)

    pmf   = pdf.(estimator.XX[:,1], estimator.XX[:,2], estimator.design, p)
    return sum( (estimator.bounds[:,2] .- estimator.bounds[:,1]) .* pmf )
end
