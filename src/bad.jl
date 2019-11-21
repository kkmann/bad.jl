module bad

import Base.show, Base.isless, Base.isequal, Base.-, Base.+, Base.*, Base.string

import Printf.@printf, Printf.@sprintf

import Distributions, Distributions.pdf, Distributions.cdf

import QuadGK.quadgk, QuadGK.gauss

using JuMP, GLPK
GLPK.jl_set_preemptive_check(false) # faster!

include("util.jl")
EarlyFutility, EarlyEfficacy = Futility(), Efficacy()
export Futility, Efficacy, CriticalValue, valid, EarlyFutility, EarlyEfficacy

include("Prior.jl")
export is_proper, condition, update, predictive_pmf, mean, expected_value

include("Beta.jl")
export Beta

include("Design.jl")
export Design, n, n1, n2, c2, power, probability, reject_null

include("optimize.jl")
export get_optimal_design

end # module
