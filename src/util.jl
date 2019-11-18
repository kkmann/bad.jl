abstract type EarlyStopping end

primitive type Futility <: EarlyStopping 8 end
Futility() = reinterpret(Futility, Int8(0))
isless(x::Int, y::Futility) = false
isless(x::Futility, y::Int) = true
isequal(x::Int, y::Futility) = false
isequal(x::Futility, y::Int) = false
-(x::Int, y::Futility) = Inf
-(x::Futility, y::Int) = -Inf


primitive type Efficacy <: EarlyStopping 8 end
Efficacy() = reinterpret(Efficacy, Int8(0))
isless(x::Int, y::Efficacy) = true
isless(x::Efficacy, y::Int) = false
isequal(x::Int, y::Efficacy) = false
isequal(x::Efficacy, y::Int) = false
-(x::Int, y::Efficacy) = -Inf
-(x::Efficacy, y::Int) = Inf

CriticalValue = Union{Int, Futility, Efficacy}

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
