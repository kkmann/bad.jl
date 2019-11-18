module bad

import Base.show, Base.isless, Base.isequal, Base.-

import Printf.@printf

import Distributions.Binomial, Distributions.Beta, Distributions.pdf, Distributions.cdf
dbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = pdf(Binomial(n, p), x)
pbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = cdf(Binomial(n, p), x)
dbeta(p::T, a::T, b::T) where{T<:Real} = pdf(Beta(a, b), p)
pbeta(p::T, a::T, b::T) where{T<:Real} = cdf(Beta(a, b), p)

import QuadGK.quadgk, QuadGK.gauss

export valid, probability

include("Prior.jl")
export Prior, condition, posterior

include("util.jl")
export Futility, Efficacy, valid

include("Design.jl")
export Design, n, n1, n2, c, power, probability, reject_null, get_x1_x2_grid

end # module
