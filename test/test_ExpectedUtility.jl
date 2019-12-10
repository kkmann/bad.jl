pnull, pmcr, palt = .2, .25, .4
α, β = .05, .2
null = PointMass(pnull)
mtoer = TypeOneErrorRate(null, pnull)
prior1 = Beta(mean = .4, sd = .1)
prior2 = .8*prior1 + .2*Beta(1, 1)
prior3 = prior2 <= .6
prior = update(prior3, 4, 10)

av_toer = TypeOneErrorRate(prior, pnull)
power = Power(prior, pmcr)

evaluate(.1*power, design)

PoS = (1 - cdf(pmcr, prior))*power
PoTOE = cdf(pnull, prior)*av_toer
ess = SampleSize(prior)

utility = (-600*PoTOE) + (300*PoS) - (1*ess)

tmp = convert(CompositeScore, utility)

minimise(utility)

utility(3, 10, 20, 13.)

minimise(ess)(3, 10, 20, 13.)

problem = Problem(
    minimise(utility),
    subject_to(mtoer, α),
    subject_to(power, β,)
)
design = optimise(problem)

plot(design)

import Gadfly
p = 0:.01:1
Gadfly.plot(x = collect(p), y = power.(design, p))

power(design, .2)
power(design, .3)
power(design, .4)
