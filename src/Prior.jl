struct Prior{T<:Real}
    pivots::Vector{T}
    weights::Vector{T}
    pdf::Vector{T}
    cdf::Vector{T}
    a::T
    b::T
    low::T
    high::T
end

# make priors iterable
Base.iterate(design::Prior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Prior) = 1

function Prior(a::Real, b::Real; low::Real = 0, high::Real = 1)
    p, ω = gauss_legendre_25(low, high)
    pdf  = dbeta.(p, convert(Float64, a), convert(Float64, b))
    z    = sum(pdf .* ω)
    pdf  = pdf ./ z # normalize
    cdf  = cumsum(pdf)
    cdf  = cdf ./ cdf[end] # normalize
    return Prior{Float64}(
        convert(Vector{Float64}, p),
        convert(Vector{Float64}, ω),
        convert(Vector{Float64}, pdf),
        convert(Vector{Float64}, cdf),
        convert(Float64, a),
        convert(Float64, b),
        convert(Float64, low),
        convert(Float64, high)
    )
end

Base.show(io::IO, prior::Prior) = @printf "Beta(a=%.2f,b=%.2f)[%.2f,%.2f]" prior.a prior.b prior.low prior.high

integrate(prior::Prior, values_on_pivots) = sum( values_on_pivots .* prior.weights )

condition(prior::Prior{T}; low::T = prior.low, high::T = prior.high) where {T<:Real} =
    Prior(prior.a, prior.b, low = max(low, prior.low), high = min(high, prior.high))

function posterior(prior::Prior{T}, x::Int, n::Int) where {T<:Real}
    !(0 <= x <= n) ? error("invalid x / n") : nothing
    Prior(prior.a + x, prior.b + n - x, low = prior.low, high = prior.high)
end
