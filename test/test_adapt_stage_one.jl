using Test

prior = Beta(mean = .35, sd = .1)
pnull = .2
pmcr  = .3
α, β  = .05, .2

design = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, pmcr, 1 - β)
    ) |>
    optimise

samplesize = SampleSize(prior)
pow        = Power(prior, pmcr)
toer       = TypeOneErrorRate(PointMass(pnull), pnull)

evaluate(samplesize, design)
evaluate(pow, design)
evaluate(toer, design)


using JuMP, GLPK

m = Model()
