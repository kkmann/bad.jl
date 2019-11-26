module bad

import Base.show, Base.isless, Base.isequal, Base.-, Base.+, Base.*, Base.string,
    Base.convert

import ProgressMeter.@showprogress

import Printf.@printf, Printf.@sprintf

import SpecialFunctions.gamma, SpecialFunctions.beta_inc

import Distributions

import QuadGK.quadgk, QuadGK.gauss

using JuMP, GLPK
GLPK.jl_set_preemptive_check(false) # faster!

import DataFrames

include("priors/Prior.jl")
export is_proper, condition, update, predictive_pmf, mean, expected_value

include("util.jl")
export valid

include("priors/Beta.jl")
export Beta

include("priors/PointMass.jl")
export PointMass



include("designs.jl")
export Design, OptimalDesign, n1, n2, c2, early_futility, early_efficacy, probability,
    reject_null, sample_space

include("power.jl")
export power



include("constraints/constraints.jl")
include("objectives/Objective.jl")

include("Problem.jl")
export Problem, OptimalDesign, optimise



include("constraints/ExpectedPowerConstraint.jl")
export minimal_expected_power

include("constraints/MaximalTypeOneErrorRateConstraint.jl")
export maximal_type_one_error_rate



include("objectives/ExpectedSampleSize.jl")
export minimise_expected_sample_size



include("estimators/Estimator.jl")
export bias, mean_squared_error

include("estimators/MaximumLikelihoodEstimator.jl")
export MaximumLikelihoodEstimator

include("estimators/PosteriorMeanEstimator.jl")
export PosteriorMeanEstimator



include("orderings/Ordering.jl")
export smaller_or_equal, strictly_smaller, larger_or_equal, strictly_larger, p_value

include("orderings/EstimatorOrdering.jl")
export EstimatorOrdering



include("confidence-intervals/ConfidenceInterval.jl")
export ConfidenceInterval

include("confidence-intervals/ClopperPearsonInterval.jl")
export ClopperPearsonInterval

end # module
