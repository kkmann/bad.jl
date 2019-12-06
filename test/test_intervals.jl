using Test

p0, p1 = .2,  .4
α, β   = .05, .2

design = Problem(
    minimise_expected_sample_size(PointMass(p1)),
    maximal_type_one_error_rate(p0, α),
    minimal_expected_power(PointMass(p1), p0, 1 - β),
) |> optimise

mle      = MaximumLikelihoodEstimator()
ordering = EstimatorOrdering(mle)

ci = ClopperPearsonInterval(ordering, design, α)

ci([0], [0])
coverage_probability(ci, [.5])
mean_width(ci, .3)



@test compatible(ci, design, p0)

prior = Beta(5, 7)
cci = PosteriorCredibleInterval(prior, design, α)

cci(0, 0)
coverage_probability(cci, .5)
mean_width(cci, .3)

@test !compatible(cci, design, p0)

jprior = JeffreysPrior(design)
jci = PosteriorCredibleInterval(jprior, design, α)

cci(0, 0)
coverage_probability(cci, .5)
mean_width(cci, .3)

@test !compatible(cci, design, p0)
