struct JeffreysPrior{T<:Real,TD<:AbstractDesign} <: Prior
    pivots::Vector{T}
    weights::Vector{T}
    design::TD
    pdf::Vector{T}
    pdf_function::Function
    low::T
    high::T
end

function JeffreysPrior(design::TD; low::Real = 0, high::Real = 1)::JeffreysPrior{Float64,TD} where
    {TD<:AbstractDesign}

    eps          = 1e-3
    low, high    = max(low, eps), min(high, 1 - eps)
    p, ω         = gauss_legendre_25(low, high)
    z            = quadgk( p -> fisher_information(p, design), low, high; atol = 1e-5 )[1]
    pdf_function = p -> ( (p < low) | (p > high) ) ? 0.0 : fisher_information(p, design) / z
    pdf          = pdf_function.(p)
    zz           = (high - low)/2 * sum(pdf .* ω)
    pdf          = pdf ./ zz # normalize
    return JeffreysPrior{Float64,TD}(
        convert(Vector{Float64}, p),
        convert(Vector{Float64}, ω),
        design,
        convert(Vector{Float64}, pdf),
        pdf_function,
        convert(Float64, low),
        convert(Float64, high)
    )
end

function pdf(p::TR, prior::JeffreysPrior{TR,TD})::TR where
    {TR<:Real,TD<:AbstractDesign}

    return prior.pdf_function(p)
end

function cdf(p::TR, prior::JeffreysPrior{TR,TD})::TR where
    {TR<:Real,TD<:AbstractDesign}

    return min(1.0, max(0.0, quadgk(prior.pdf_function, 0, p; atol = 1e-5)[1]))
end

function expectation(f::Function, prior::JeffreysPrior{TR,TD})::TR where
    {TR<:Real,TD<:AbstractDesign}

    return (prior.high - prior.low)/2 * sum( prior.pdf .* f.(prior.pivots) .* prior.weights )
end

function mean(prior::JeffreysPrior{TR,TD})::TR where
    {TR<:Real,TD<:AbstractDesign}

    return expectation(p -> p, prior)
end

function condition(prior::JeffreysPrior{TR,TD}; low::TR = prior.low, high::TR = prior.high)::JeffreysPrior{TR,TD} where
    {TR<:Real,TD<:AbstractDesign}

    return JeffreysPrior{TR,TD}(TD; low = max(low, prior.low), high = min(high, prior.high))
end

function update(prior::JeffreysPrior{TR,TD}, x::TI, n::TI)::GenericDistribution{TR} where
    {TR<:Real,TI<:Integer,TD<:AbstractDesign}

    return GenericDistribution( p -> pdf(p, prior) * dbinom(x, n, p) )
end

function string(prior::JeffreysPrior)
    eps          = 1e-3
    low, high    = max(prior.low, eps), min(prior.high, 1 - eps)
    if (prior.low > low) | (prior.high < high)
        @sprintf "JeffreysPrior|[%.2f,%.2f](%s)" prior.low prior.high prior.design
    else
        @sprintf "JeffreysPrior(%s)" prior.design
    end
end
