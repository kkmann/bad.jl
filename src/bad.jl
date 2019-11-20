module bad

import Base.show, Base.isless, Base.isequal, Base.-

import Printf.@printf

import Distributions, Distributions.pdf, Distributions.cdf

import QuadGK.quadgk, QuadGK.gauss

using JuMP, GLPK
GLPK.jl_set_preemptive_check(false) # faster!

include("util.jl")
EarlyFutility, EarlyEfficacy = Futility(), Efficacy()
export Futility, Efficacy, CriticalValue, valid, EarlyFutility, EarlyEfficacy

include("Beta.jl")
export Beta, condition, update, predictive_pmf

include("Design.jl")
export Design, n, n1, n2, c, power, probability, reject_null

include("optimize.jl")
export get_optimal_design

end # module
