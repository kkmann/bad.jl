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

@time design = optimise(problem; verbosity = 3)
power(design), mtoer(design), ess(design)

as_table(design)
design.n2


mtoer(design, 5; partial_stage_two = (2, 5))


α_new = mtoer(design; partial_stage_one = (0, 0) )
β_new = 1 - power(design; partial_stage_one = (0, 0) )

obs = 5, 10
α_new = mtoer(design; partial_stage_one = obs )
β_new = 1 - power(design; partial_stage_one = obs )

adapt(design, prior, obs)
