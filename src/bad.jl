module bad

import Base.show, Base.isless, Base.isequal, Base.-

import Distributions.Binomial, Distributions.Beta, Distributions.pdf, Distributions.cdf
dbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = pdf(Binomial(n, p), x)
pbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = cdf(Binomial(n, p), x)
dbeta(p::T, a::T, b::T) where{T<:Real} = pdf(Beta(a, b), p)
pbeta(p::T, a::T, b::T) where{T<:Real} = cdf(Beta(a, b), p)

import QuadGK.quadgk

export valid, probability

include("util.jl")
export Futility, Efficacy, valid

include("Prior.jl")
export Prior, BetaPrior, condition


include("Design.jl")
export Design, n, n1, n2, c, power

end # module
