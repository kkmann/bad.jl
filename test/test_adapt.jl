using Test, bad; import Plots

# construct the usual prior
pnull  = .2 # response under TAU
pmcr   = pnull + .1 # minimal clinically relevant
prior1 = Beta(mean = .35, sd = .1)  # start with an subjective prior
prior2 = update(prior1, 4, 10)      # update with phase I data
prior3 = .8*prior2 + .2*Beta(1, 1)  # robustify
prior  = prior3 <= min(2*pmcr, 1.0) # restrict to plausible range

mtoer  = TypeOneErrorRate(prior | pnull)
power = Power(prior >= pmcr)
ess   = SampleSize(prior)

# first consider the 'standard design'
α, β = .05, .2
problem = Problem(
    minimise(ess),
    subject_to(mtoer, α),
    subject_to(power, β)
)

design = optimise(problem; verbosity = 0)
power(design), mtoer(design), ess(design)

# test invariance
obs = 0, 0
α_new = mtoer(design; partial_stage_one = obs )
β_new = 1 - power(design; partial_stage_one = obs )

adesign = adapt(design, prior, obs)
@test all(adesign.n2 .== design.n2)
@test all(adesign.c2 .== design.c2)
@test power(design) ≈ power(adesign; partial_stage_one = obs)
@test mtoer(design) ≈ mtoer(adesign; partial_stage_one = obs)
@test ess(design)   ≈ ess(adesign; partial_stage_one = obs)

obs   = 5, 10
α_new = mtoer(design; partial_stage_one = obs )
β_new = 1 - power(design; partial_stage_one = obs )
adesign = adapt(design, prior, obs)

@test power(adesign; partial_stage_one = obs) >= 1 - β_new
@test mtoer(adesign; partial_stage_one = obs) <= α_new
@test ess(adesign) >= ess(design)
@test ess(adesign; partial_stage_one = obs) >= ess(design; partial_stage_one = obs)


prior2 = update(prior, 6, 10)
adesign = adapt(design, prior2, obs)
@test Power(prior2 >= pmcr)(adesign; partial_stage_one = obs) >= 1 - β_new
@test TypeOneErrorRate(prior2 | pnull)(adesign; partial_stage_one = obs) <= α_new
