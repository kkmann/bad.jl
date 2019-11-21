struct Beta{T<:Real} <: Prior
    pivots::Vector{T}
    weights::Vector{T}
    pdf::Vector{T}
    a::T
    b::T
    low::T
    high::T
end

function Beta(a::Real, b::Real; low::Real = 0, high::Real = 1)
    p, ω = gauss_legendre_25(low, high)
    pdf  = dbeta.(p, convert(Float64, a), convert(Float64, b))
    z    = (high - low)/2 * sum(pdf .* ω)
    pdf  = pdf ./ z # normalize
    return Beta{Float64}(
        convert(Vector{Float64}, p),
        convert(Vector{Float64}, ω),
        convert(Vector{Float64}, pdf),
        convert(Float64, a),
        convert(Float64, b),
        convert(Float64, low),
        convert(Float64, high)
    )
end
# convenience function to fit a and b
function Beta(;mean::Real = .5, sd::Real = sqrt(1/12))
    sd >= .5 ? error("sd must be strictly smaller than .5") : nothing
    # solve for a, b
    a = (-mean^3 + mean^2 - mean*sd^2) / sd^2
    b = (mean - 1) * (mean^2 - mean + sd^2) / sd^2
    Beta(a, b; low = 0, high = 1)
end

function Base.show(io::IO, prior::Beta)
    if (prior.low != 0) | (prior.high != 1)
        @printf "Beta|[%.2f,%.2f](a=%.2f,b=%.2f)" prior.low prior.high prior.a prior.b
    else
        @printf "Beta(a=%.2f,b=%.2f)" prior.a prior.b
    end
end

condition(prior::Beta{T}; low::T = prior.low, high::T = prior.high) where {T<:Real} =
    Beta(prior.a, prior.b, low = max(low, prior.low), high = min(high, prior.high))

function update(prior::Beta{T}, x::Int, n::Int) where {T<:Real}
    !(0 <= x <= n) ? error("invalid x / n") : nothing
    Beta(prior.a + x, prior.b + n - x, low = prior.low, high = prior.high)
end

integrate(prior::Beta, values_on_pivots) = (prior.high - prior.low) / 2 * sum( prior.pdf .* values_on_pivots .* prior.weights )

mean(prior::Beta) = integrate(prior, prior.pivots)

function predictive_pmf(x, n, prior::Beta{T}) where {T<:Real}
    integrate(prior, dbinom.(x, n, prior.pivots))
end
