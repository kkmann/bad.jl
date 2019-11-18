module bad

import Base.show
import Distributions.Binomial, Distributions.pdf
dbinom(x::TI, n::TI, p::TR) where{TI<:Integer, TR<:Real} = pdf(Binomial(n, p), x)

export valid, probability

include("util.jl")
export Futility, Efficacy, valid

include("Design.jl")
export Design, n, n1, n2, c

end # module
