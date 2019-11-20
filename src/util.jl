dbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = pdf(Distributions.Binomial(n, p), x)
pbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = cdf(Distributions.Binomial(n, p), x)
dbeta(p::T, a::T, b::T) where{T<:Real} = pdf(Distributions.Beta(a, b), p)
pbeta(p::T, a::T, b::T) where{T<:Real} = cdf(Distributions.Beta(a, b), p)

abstract type EarlyStopping end
# make iterable:
Base.iterate(design::EarlyStopping, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::EarlyStopping) = 1

primitive type Futility <: EarlyStopping 8 end
Futility() = reinterpret(Futility, Int8(0))
isless(x::Int, y::Futility) = true
isless(x::Futility, y::Int) = false
isequal(x::Int, y::Futility) = false
isequal(x::Futility, y::Int) = false
-(x::Int, y::Futility) = -Inf
-(x::Futility, y::Int) = Inf


primitive type Efficacy <: EarlyStopping 8 end
Efficacy() = reinterpret(Efficacy, Int8(0))
isless(x::Int, y::Efficacy) = false
isless(x::Efficacy, y::Int) = true
isequal(x::Int, y::Efficacy) = false
isequal(x::Efficacy, y::Int) = false
-(x::Int, y::Efficacy) = Inf
-(x::Efficacy, y::Int) = -Inf

CriticalValue = Union{Int, Futility, Efficacy}

early_stop = c2 -> isa(c2, EarlyStopping)

valid(p::T) where {T<:Real} = 0 <= p <= 1

gauss_legendre = function(low, high, order::Integer)
    pivots, weights = gauss(order)
    a, b = (high - low)/2, (high + low)/2
    scaled_pivots  = a.*pivots .+ b
    scaled_weights = a.*weights
    return(scaled_pivots, scaled_weights)
end

# precompute
gl_25_pivots, gl_25_weights = gauss_legendre(-1, 1, 25)
gauss_legendre_25 = function(low, high)
    pivots, weights = gl_25_pivots, gl_25_weights
    a, b = (high - low)/2, (high + low)/2
    scaled_pivots  = a.*pivots .+ b
    scaled_weights = a.*weights
    return(scaled_pivots, scaled_weights)
end

function guess_nmax(prior, p0, mrv, α, β; multiple = 2)
    cprior  = condition(prior, low = mrv)
    # heuristic for nmax: 2* sample size of fixed z test powered for prior mean
    p1      = mean(cprior)
    z_1_α   = Distributions.quantile(Distributions.Normal(), 1 - α)
    z_1_β   = Distributions.quantile(Distributions.Normal(), 1 - β)
    napprox = p1*(1 - p1)*( (z_1_α + z_1_β) / (p1 - p0) )^2
    return Int(ceil(multiple * napprox))
end

function to_numeric(c2)
    c2 == EarlyFutility ? (return Inf) : nothing
    c2 == EarlyEfficacy ? (return -Inf) : nothing
    return convert(Float64, c2)
end
