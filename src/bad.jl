module bad

import Base.show, Base.isless, Base.isequal, Base.-
import Distributions.Binomial, Distributions.pdf, Distributions.cdf
dbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = pdf(Binomial(n, p), x)
pbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = cdf(Binomial(n, p), x)

export valid, probability

include("util.jl")
export Futility, Efficacy, valid

include("Design.jl")
export Design, n, n1, n2, c, power

end # module
