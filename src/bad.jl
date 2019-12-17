module bad

import Base.show, Base.isless, Base.isequal, Base.-, Base.+, Base.*, Base.string,
    Base.convert, Base.<=, Base.>=, Base.|, Base.length, Base.size

import Printf.@printf, Printf.@sprintf

import ProgressMeter

import SpecialFunctions.gamma, SpecialFunctions.beta_inc

import Distributions, Distributions.cdf, Distributions.pdf

import QuadGK.quadgk, QuadGK.gauss

import Roots

using JuMP, GLPK, Ipopt
GLPK.jl_set_preemptive_check(false) # faster!

import DataFrames



include("priors/Prior.jl")
export is_proper, condition, update, pdf, cdf, mean, expectation

include("util.jl")
export valid



include("designs.jl")
export Design, OptimalDesign,
    n1, n2, n, c2,
    early_futility, early_efficacy, continuation_region, futility_region,
    efficacy_region, early_stop_region,
    as_table, reject, sample_space



include("priors/JeffreysPrior.jl")
export JeffreysPrior

include("priors/GenericDistribution.jl")
export GenericDistribution

include("priors/Beta.jl")
export Beta, BetaMixture

include("priors/PointMass.jl")
export PointMass



include("pmf.jl")
export pmf, pmf_x1, pmf_x2_given_x1, cdf, cdf_x1, cdf_x2_given_x1



include("estimators/Estimator.jl")
export estimate, bias, mean_squared_error, mean_absolute_error

include("estimators/MaximumLikelihoodEstimator.jl")
export MaximumLikelihoodEstimator

include("estimators/PosteriorMean.jl")
export PosteriorMean, PosteriorMeanPrecalculated

include("estimators/RaoBlackwellEstimator.jl")
export RaoBlackwellEstimator

include("estimators/CompatibleMLE.jl")
export CompatibleMLE



include("orderings/Ordering.jl")
export EstimatorOrdering, more_extreme, p_value,
    compatible, mlecompatible

include("orderings/PValue.jl")
export PValue, evaluate



include("interval-estimators/IntervalEstimator.jl")
export IntervalEstimator, get_bounds, coverage_probability, mean_width,
    compatible

include("interval-estimators/ClopperPearsonInterval.jl")
export ClopperPearsonInterval

include("interval-estimators/PosteriorCredibleInterval.jl")
export PosteriorCredibleInterval



include("Score.jl")
export update!,
    SampleSize, Power, TypeOneErrorRate,
    CompositeScore


include("Problem.jl")
export Problem, optimise, adapt


include("objectives/Objective.jl")
export ScoreObjective, minimise, maximise

include("objectives/MiniMaxSampleSize.jl")
export MiniMaxSampleSize

include("constraints.jl")
export PowerConstraint, TypeOneErrorRateConstraint, subject_to






end # module
