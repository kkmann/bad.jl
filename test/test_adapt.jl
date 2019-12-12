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

aproblem = AdaptationProblem(design, 6, 15)

tmp, ind, n1_selected = bad.build_model(aproblem)

bad.optimise!(tmp, 3, 60)


aproblem2 = AdaptationProblem(design, 0, 0)

tmp, ind, n1_selected = bad.build_model(aproblem2)

bad.optimise!(tmp, 3, 60)


plot(design)

power(design)

power(design; x1partial = 6, n1partial = 15)
mtoer(design; x1partial = 6, n1partial = 15)

6/15

@test begin
    map(
        x1p -> power(design; x1partial = x1p, n1partial = 16),
        0:15
    ) |>
    x -> diff(x) |>
    x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
end
@test begin
    map(
        x1p -> toer(design; x1partial = x1p, n1partial = 16),
        0:15
    ) |>
    x -> diff(x) |>
    x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
end


xx1, nn2, cc2 = bad.extract_solution(aproblem, ind, n1_selected)

@time design = optimise(aproblem; verbosity = 3)
power(design), mtoer(design), ess(design)

bad.grid(problem)
bad.grid(aproblem)
