using Test, bad

import Plots



include("test_util.jl")

include("test_Beta.jl")
include("test_BetaMixture.jl")



# define a few standards
pnull, pmcr, palt = .2, .25, .4
α, β = .05, .2
null = PointMass(pnull)
mtoer = TypeOneErrorRate(null, pnull)
prior1 = Beta(mean = .4, sd = .1)
prior2 = .8*prior1 + .2*Beta(1, 1)
prior3 = prior2 <= .6
prior = update(prior3, 4, 10)
p = 0:.01:1
pdfs   = [pdf.(p, prior) for prior in (prior1, prior2, prior3, prior)] |> x -> hcat(x...)
cdfs   = [cdf.(p, prior) for prior in (prior1, prior2, prior3, prior)] |> x -> hcat(x...)
Plots.plot(p, pdfs)
Plots.plot(p, cdfs)



include("test_Simon.jl")



power = Power(prior, pmcr)
mean(prior), mean(power.prior)
problem = Problem(
    minimise(SampleSize(prior)),
    subject_to(mtoer, α),
    subject_to(power, β,)
)
design = optimise(problem)
plot(design)



include("test_pmf.jl")

include("test_scores.jl")

include("test_JeffreysPrior.jl")



include("test_PosteriorMean.jl")

include("test_RaoBlackwellEstimator.jl")

include("test_CompatibleMLE.jl")



power = Power(prior, pmcr)
mean(prior), mean(power.prior)
problem = Problem(
    minimise(SampleSize(prior)),
    subject_to(mtoer, α),
    subject_to(power, β,)
)
design = optimise(problem)
plot(design)



include("test_intervals.jl")



include("test_MiniMaxSampleSize.jl")
