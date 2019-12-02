module bad

import Base.show, Base.isless, Base.isequal, Base.-, Base.+, Base.*, Base.string,
    Base.convert, Base.<=, Base.>=, Base.|

import ProgressMeter

import TickTock.tick, TickTock.tok

import Printf.@printf, Printf.@sprintf

import SpecialFunctions.gamma, SpecialFunctions.beta_inc

import Distributions, Distributions.cdf, Distributions.pdf

import QuadGK.quadgk, QuadGK.gauss

using JuMP, GLPK
GLPK.jl_set_preemptive_check(false) # faster!

import DataFrames, Gadfly

include("priors/Prior.jl")
export is_proper, condition, update, pdf, cdf, mean, expectation

include("util.jl")
export valid


include("designs.jl")
export Design, OptimalDesign, n1, n2, n, c2, early_futility, early_efficacy, as_table,
    reject_null, sample_space, plot, expected_sample_size

include("power.jl")
export power

include("priors/JeffreysPrior.jl")
export JeffreysPrior

include("priors/GenericDistribution.jl")
export GenericDistribution

include("priors/Beta.jl")
export Beta, BetaMixture

include("priors/PointMass.jl")
export PointMass


include("constraints/constraints.jl")
export update!

include("objectives/Objective.jl")

include("Problem.jl")
export Problem, OptimalDesign, optimise


include("constraints/no-constraints.jl")
export NoPowerConstraint, NoTypeOneErrorRateConstraint

include("constraints/ExpectedPowerConstraint.jl")
export minimal_expected_power

include("constraints/MaximalTypeOneErrorRateConstraint.jl")
export maximal_type_one_error_rate



include("objectives/ExpectedSampleSize.jl")
export minimise_expected_sample_size

include("objectives/MiniMaxSampleSize.jl")
export MiniMaxSampleSize

include("objectives/ExpectedUtility.jl")
export ExpectedUtility


include("adapt.jl")
export adapt



include("estimators/Estimator.jl")
export estimate, bias, mean_squared_error, mean_absolute_error

include("estimators/MaximumLikelihoodEstimator.jl")
export MaximumLikelihoodEstimator

include("estimators/PosteriorMean.jl")
export PosteriorMean



include("orderings/Ordering.jl")
export smaller_or_equal, strictly_smaller, larger_or_equal, strictly_larger, p_value

include("orderings/EstimatorOrdering.jl")
export EstimatorOrdering



include("confidence-intervals/ConfidenceInterval.jl")
export ConfidenceInterval

include("confidence-intervals/ClopperPearsonInterval.jl")
export ClopperPearsonInterval

end # module
