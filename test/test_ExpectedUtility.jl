using Test, bad

# construct the usual prior
pnull  = .2 # response under TAU
pmcr   = pnull + .1 # minimal clinically relevant
prior1 = Beta(mean = .35, sd = .1)  # start with an subjective prior
prior2 = update(prior1, 4, 10)      # update with phase I data
prior3 = .8*prior2 + .2*Beta(1, 1)  # robustify
prior  = prior3 <= min(2*pmcr, 1.0) # restrict to plausible range

p = 0:.01:1

# probability of succes is simply power * Pr[ p >= pmcr ]
pos   =  (1 - cdf(pmcr, prior)) * Power(prior >= pmcr)
# probability of an type one error is avtoer * Pr [p <= pnull ]
potoe = cdf(pnull, prior) * Power(prior <= pnull)
# finally, expected sample size is the usual
ess   = SampleSize(prior)
# overall utility is then goverend by two factors, both on the scale of the
# average per-patient costs in phase II
var"cost failed phase III" = 12000.
var"risk weighted profit of successful phaseIII" = 200.
utility = -var"cost failed phase III"*potoe +
    var"risk weighted profit of successful phaseIII"*pos -
    ess

# first consider the 'standard design'
α, β = .05, .2
problem = Problem(
    minimise(ess),
    Power(prior  | pnull) <= α,
    Power(prior >= pmcr)  >= 1 - β
)

design = optimise(problem; verbosity = 0)
Power(prior >= pmcr)(design), Power(prior | pnull)(design), ess(design), utility(design)

# now, lets relax the power and type one error rate constraints and maximise
# utility instead! Note that we need to manually increase the marginal
# feasible space since the heuristics based on the given constraints are
# now favouring small designs
uproblem = Problem(
    maximise(utility),
    Power(prior  | pnull) <= 4*α,
    Power(prior >= pmcr)  >= 1 - 2*β,
    n1values = 5:35,
    nmax     = 100
)
udesign = optimise(uproblem; verbosity = 0)
Power(prior >= pmcr)(udesign), Power(prior | pnull)(udesign), ess(udesign), utility(udesign)

@test utility(udesign) > utility(design)
