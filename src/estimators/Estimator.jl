abstract type Estimator end

Base.iterate(estimator::Estimator, state = 0) = state > 0 ? nothing : (estimator, state + 1)
Base.length(estimator::Estimator) = 1

Base.show(io::IO, estimator::Estimator) = print(io, string(estimator))
Base.show(io::IO, ::MIME"application/prs.juno.inline", estimator::Estimator) = print(io, string(estimator))

function bias(p::Real, estimator::Estimator, design::AbstractDesign)
    0 <= p <= 1 ? nothing : error("p must be between 0 and 1")
    space = sample_space(design)
    return (estimator.(space[:,1], space[:,2], design) .- p) .* probability.(space[:,1], space[:,2], design, p) |> sum
end

function mean_squared_error(p::Real, estimator::Estimator, design::AbstractDesign)
    0 <= p <= 1 ? nothing : error("p must be between 0 and 1")
    space = sample_space(design)
    return (estimator.(space[:,1], space[:,2], design) .- p).^2 .* probability.(space[:,1], space[:,2], design, p) |> sum
end
