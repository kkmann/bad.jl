# function dbinom(x::TI, n::TI, p::TR) where {TI<:Integer,TR<:Real}
#     return (0 <= x <= n) ? gamma(n + 1)/gamma(x + 1)/gamma(n - x + 1)*(1 - p)^(n - x)*p^x : 0.0
# end
# function dbinom(x::TI, n::TI, prior::TP) where {TI<:Integer,TP<:Prior}
#     expectation(p -> dbinom(x, n, p), prior)
# end
#
# function pbinom(x::Real, n::Real, p::Real)
#     n < 0  ? (return NaN) : nothing
#     x < 0  ? (return 0.0) : nothing
#     x >= n ? (return 1.0) : nothing
#     return beta_inc(n - x, x + 1, 1 - p, p)[1]
# end
# pbinom(x::Real, n::Real, p::Prior) = expectation(p -> pbinom(x, n, p), p)

function dbeta(p::Real, a::Real, b::Real)
    any((a, b) .<= 0) ? (return NaN) : nothing
    !(0 <= p <= 1) ? (return 0.0) : nothing
    return gamma(a + b)/gamma(a)/gamma(b)*(1 - p)^(b - 1)*p^(a - 1)
end
function pbeta(p::Real, a::Real, b::Real)
    any((a, b) .<= 0) ? (return NaN) : nothing
    p <= 0  ? (return 0.0) : nothing
    p >= 1  ? (return 1.0) : nothing
    return beta_inc(a, b, p, 1 - p)[1]
end

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

function one_stage_sample_size(p0::Real, α::Real, p1::Real, β::Real)
    z_1_α   = Distributions.quantile(Distributions.Normal(), 1 - α)
    z_1_β   = Distributions.quantile(Distributions.Normal(), 1 - β)
    napprox = p1*(1 - p1)*( (z_1_α + z_1_β) / (p1 - p0) )^2
    return Int(ceil(napprox))
end
