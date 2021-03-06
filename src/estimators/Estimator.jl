abstract type Estimator end

Base.iterate(estimator::Estimator, state = 0) = state > 0 ? nothing : (estimator, state + 1)
Base.length(estimator::Estimator) = 1

Base.show(io::IO, estimator::Estimator) = print(io, string(estimator))
Base.show(io::IO, ::MIME"application/prs.juno.inline", estimator::Estimator) = print(io, string(estimator))

(::Estimator)(x1::Integer, x2::Integer, design::AbstractDesign) = error("not implemented")
function estimate(estimator::TE, x1::TI, x2::TI, design::TD) where {TE<:Estimator,TI<:Integer,TD<:AbstractDesign}
    return estimator(x1, x2, design)
end


string(estimator::Estimator) = error("not implemented")

function bias(p::Real, estimator::Estimator, design::AbstractDesign)
    0 <= p <= 1 ? nothing : error("p must be between 0 and 1")
    XX = sample_space(design)
    x1, x2 = XX[:,1], XX[:,2]
    return (estimator.(x1, x2, design) .- p) .*
        pmf.(x2, n2.(design, x1), p) .*
        pmf.(x1, n1(design), p) |>
        sum
end

function mean_squared_error(p::Real, estimator::Estimator, design::AbstractDesign)
    0 <= p <= 1 ? nothing : error("p must be between 0 and 1")
    XX = sample_space(design)
    x1, x2 = XX[:,1], XX[:,2]
    return (estimator.(x1, x2, design) .- p).^2 .*
        pmf.(x2, n2.(design, x1), p) .*
        pmf.(x1, n1(design), p) |>
        sum
end

function mean_absolute_error(p::Real, estimator::Estimator, design::AbstractDesign)
    0 <= p <= 1 ? nothing : error("p must be between 0 and 1")
    XX = sample_space(design)
    x1, x2 = XX[:,1], XX[:,2]
    return abs.(estimator.(x1, x2, design) .- p) .*
        pmf.(x2, n2.(design, x1), p) .*
        pmf.(x1, n1(design), p) |>
        sum
end
