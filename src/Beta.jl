struct Beta{T<:Real}
    pivots::Vector{T}
    weights::Vector{T}
    pdf::Vector{T}
    cdf::Vector{T}
    a::T
    b::T
    low::T
    high::T
end

# make betas iterable
Base.iterate(design::Beta, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Beta) = 1

function Beta(a::Real, b::Real; low::Real = 0, high::Real = 1)
    p, ω = gauss_legendre_25(low, high)
    pdf  = dbeta.(p, convert(Float64, a), convert(Float64, b))
    pdf  = pdf ./ sum(pdf .* ω) # normalize
    cdf  = cumsum(pdf)
    cdf  = cdf ./ cdf[end] # normalize
    return Beta{Float64}(
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
# convenience function to fit a and b
function Beta(;mean::Real = .5, sd::Real = sqrt(1/12))
    sd >= .5 ? error("sd must be strictly smaller than .5") : nothing
    # solve for a, b
    a = (-mean^3 + mean^2 - mean*sd^2) / sd^2
    b = (mean - 1) * (mean^2 - mean + sd^2) / sd^2
    Beta(a, b; low = 0, high = 1)
end

function Base.show(io::IO, beta::Beta)
    if (beta.low != 0) | (beta.high != 1)
        @printf "Beta|[%.2f,%.2f](a=%.2f,b=%.2f)" beta.low beta.high beta.a beta.b
    else
        @printf "Beta(a=%.2f,b=%.2f)" beta.a beta.b
    end
end

condition(beta::Beta{T}; low::T = beta.low, high::T = beta.high) where {T<:Real} =
    Beta(beta.a, beta.b, low = max(low, beta.low), high = min(high, beta.high))

function update(beta::Beta{T}, x::Int, n::Int) where {T<:Real}
    !(0 <= x <= n) ? error("invalid x / n") : nothing
    Beta(beta.a + x, beta.b + n - x, low = beta.low, high = beta.high)
end

integrate(beta::Beta, values_on_pivots) = (beta.high - beta.low) / 2 * sum( beta.pdf .* values_on_pivots .* beta.weights )

function predictive_pmf(x, n, beta::Beta{T}) where {T<:Real}
    integrate(beta, dbinom.(x, n, beta.pivots))
end
