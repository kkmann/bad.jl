using Test, bad; import Plots



pnull = .3
pmcr  = .4
p     = Beta(6, 5)
α, β  = .05, .2
problem = Problem(
    minimise(SampleSize(p)),
    subject_to(TypeOneErrorRate(p | pnull), α),
    subject_to(Power(p >= pmcr), β)
)
design = optimise(problem; verbosity = 0)



mle      = MaximumLikelihoodEstimator()
ordering = EstimatorOrdering(mle)
ci = ClopperPearsonInterval(ordering, design, α)
ci([0], [0])
coverage_probability(ci, [.5])
mean_width(ci, .3)
@test compatible(ci, design, pnull)["compatible"]




pci = PosteriorCredibleInterval(Beta(3, 7), design, α)
pci(0, 0)
coverage_probability(pci, .5)
mean_width(pci, .3)
@test !compatible(pci, design, pnull)["compatible"]




jprior = JeffreysPrior(design)
jci = PosteriorCredibleInterval(jprior, design, α)
jci(0, 0)
coverage_probability(jci, .5)
mean_width(jci, .3)
@test compatible(jci, design, pnull)["compatible"]
