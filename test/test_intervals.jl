using Test

mle      = MaximumLikelihoodEstimator()
ordering = EstimatorOrdering(mle)

ci = ClopperPearsonInterval(ordering, design, α)

ci([0], [0])
coverage_probability(ci, [.5])
mean_width(ci, .3)



@test compatible(ci, design, pnull)["compatible"]

prior = Beta(5, 7)
cci = PosteriorCredibleInterval(prior, design, α)

cci(0, 0)
coverage_probability(cci, .5)
mean_width(cci, .3)

@test !compatible(cci, design, p0)["compatible"]

jprior = JeffreysPrior(design)
jci = PosteriorCredibleInterval(jprior, design, α)

cci(0, 0)
coverage_probability(cci, .5)
mean_width(cci, .3)

@test !compatible(cci, design, p0)["compatible"]
