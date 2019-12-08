struct GenericDistribution{T<:Real} <: Prior
    pivots::Vector{T}
    weights::Vector{T}
    pdf::Vector{T}
    pdf_function::Function
    low::T
    high::T
end

function GenericDistribution(f::Function; low::Real = 0, high::Real = 1)::GenericDistribution{Float64}

    eps          = 1e-3
    low, high    = max(low, eps), min(high, 1 - eps)
    p, ω         = gauss_legendre_25(low, high)
    z            = quadgk( p -> f(p), low, high )[1]
    pdf_function = p -> ( (p < low) | (p > high) ) ? 0.0 : f(p) / z
    pdf          = pdf_function.(p)
    zz           = (high - low)/2 * sum(pdf .* ω)
    pdf          = pdf ./ zz # normalize
    return GenericDistribution{Float64}(
        convert(Vector{Float64}, p),
        convert(Vector{Float64}, ω),
        convert(Vector{Float64}, pdf),
        pdf_function,
        convert(Float64, low),
        convert(Float64, high)
    )
end

function pdf(p::TR, prior::GenericDistribution{TR})::TR where {TR<:Real}

    return prior.pdf_function(p)
end

function cdf(p::TR, prior::GenericDistribution{TR})::TR where {TR<:Real}

    return max(0.0, min(1.0, quadgk(prior.pdf_function, 0, p)[1]))
end

function expectation(f::Function, prior::GenericDistribution{TR})::TR where
    {TR<:Real}

    return (prior.high - prior.low)/2 * sum( prior.pdf .* f.(prior.pivots) .* prior.weights )
end

function mean(prior::GenericDistribution{TR})::TR where
    {TR<:Real}

    return expectation(p -> p, prior)
end

function condition(prior::GenericDistribution{TR}; low::TR = prior.low, high::TR = prior.high)::GenericDistribution{TR} where
    {TR<:Real}

    return GenericDistribution{TR}(TD; low = max(low, prior.low), high = min(high, prior.high))
end

function update(prior::GenericDistribution{TR}, x::TI, n::TI)::GenericDistribution{TR} where
    {TR<:Real,TI<:Integer}

    return GenericDistribution{TR}( p -> pdf(prior, p) * dbinom(x, n, p); low = prior.low, high = prior.high )
end

function string(prior::GenericDistribution)
    eps          = 1e-3
    low, high    = max(prior.low, eps), min(prior.high, 1 - eps)
    if (prior.low > low) | (prior.high < high)
        @sprintf "GenericDistribution|[%.2f,%.2f](%s)" prior.low prior.high prior.design
    else
        "GenericDistribution|[0, 1]"
    end
end
