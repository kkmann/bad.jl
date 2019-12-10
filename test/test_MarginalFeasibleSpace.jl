Pkg.activate("."); using Revise
using bad, Test

prior = Beta(mean = .35, sd = .1)
pnull = .2
pmcr  = .3
α, β  = .05, .2

pow   = Power(prior, pmcr)
mtoer = TypeOneErrorRate(PointMass(pnull), pnull)

problem = Problem(
    minimise(SampleSize(prior)),
    subject_to(mtoer, α, [α/5, 5*α]),
    subject_to(pow, β, [.5, .99])
)

size(problem)

design = optimise(problem; verbosity = 0, timelimit = 60)

SampleSize(prior)(design)

design.info
