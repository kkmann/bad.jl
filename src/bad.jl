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

include("AbstractDesign.jl")
export n, n1, n2, c2, early_futility, early_efficacy, power, probability, reject_null

include("Design.jl")
export Design

include("DesignIPModel.jl")
export DesignIPModel

include("OptimalDesign.jl")
export OptimalDesign

include("optimise.jl")
export optimise

include("ExpectedPowerConstraint.jl")
export minimal_expected_power

include("ExpectedSampleSize.jl")
export minimize_expected_sample_size

end # module
