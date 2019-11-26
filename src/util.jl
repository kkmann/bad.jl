dbinom(x::Int, n::Int, p::T) where {T<:Real} = pdf(Distributions.Binomial(n, p), x)
dbinom(x::Int, n::Int, p::T) where {T<:Prior} = expected_value(p -> dbinom(x, n, p), p)
pbinom(x::Int, n::Int, p::T) where {T<:Real} = cdf(Distributions.Binomial(n, p), x)
pbinom(x::Int, n::Int, p::T) where {T<:Prior} = expected_value(p -> pbinom(x, n, p), p)
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
    scaled_pivots  = a .* pivots .+ b
    scaled_weights = a .* weights
    return(scaled_pivots, scaled_weights)
end

# precompute
gl_25_pivots, gl_25_weights = gauss_legendre(-1, 1, 25)
gauss_legendre_25 = function(low, high)
    pivots, weights = gl_25_pivots, gl_25_weights
    a, b = (high - low)/2, (high + low)/2
    scaled_pivots  = a .* pivots .+ b
    scaled_weights = a .* weights
    return(scaled_pivots, scaled_weights)
end

function to_numeric(c2)
    c2 == EarlyFutility ? (return Inf) : nothing
    c2 == EarlyEfficacy ? (return -Inf) : nothing
    return convert(Float64, c2)
end


power(n::Int, c::Futility, p::T) where {T<:Union{Real,Prior}} = 0.0
power(n::Int, c::Efficacy, p::T) where {T<:Union{Real,Prior}} = 1.0
power(n::Int, c::Int, p::T) where {T<:Real} = 1 - pbinom(c, n, p)
power(n::Int, c::Int, prior::Prior) = expected_value(p -> power(n, c, p), prior)

power(x1::Int, n1::Int, n2::Int, c2::Futility, prior::Prior) = 0.0
power(x1::Int, n1::Int, n2::Int, c2::Efficacy, prior::Prior) = 1.0
power(x1::Int, n1::Int, n2::Int, c2::Int, prior::Prior) = expected_value(p -> power(n2, c2, p), update(prior, x1, n1))

function power(p::T, xx1::Int, nn1::Int, xx2::Int, nn2::Int, n1::Int, n2::Vector{Int}, c2::Vector{Union{Int,CriticalValue}}) where {T<:Union{Real,Prior}}
    any(length.((n2, c2)) .!= n1 + 1) ? error("n2/c2 must be of length n1 + 1") : nothing
    if nn2 >= n2[xx1 + 1] # everything done
        return reject_null(xx2, c2[xx1 + 1]) ? 1.0 : 0.0
    else
        if nn2 > 0 # in stage two, nn1 = n1
            return power(n2[xx1 + 1] - nn2, c2[xx1 + 1] - xx2, p)
        else # in stage one
            x1_delta = 0:(n1 - nn1)
            x1       = xx1 .+ x1_delta
            return sum( dbinom.(x1_delta, n1 - nn1, p) .* power.(n2[x1 .+ 1], c2[x1 .+ 1], p) )
        end
    end
end

function power(p::T, x::Vector{Bool}, n1::Int, n2::Vector{Int}, c2::Vector{Union{Int,CriticalValue}}) where {T<:Union{Real,Prior}}
    nn1 = min(length(x), n1)
    xx1 = sum(x[1:nn1])
    nn2 = min(length(x) - nn1, n2[xx1 + 1])
    xx2 = sum(x[(nn1 + 1):nn2])
    return power(p, xx1, nn1, xx2, nn2, n1, n2, c2)
end


function guess_nmax(p0::Real, α::Real, p1::Real, β::Real; multiple = 2)
    z_1_α   = Distributions.quantile(Distributions.Normal(), 1 - α)
    z_1_β   = Distributions.quantile(Distributions.Normal(), 1 - β)
    napprox = p1*(1 - p1)*( (z_1_α + z_1_β) / (p1 - p0) )^2
    return Int(ceil(multiple * napprox))
end
