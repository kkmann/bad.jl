prior = Prior(mean = .4, sd = .1)

@test bad.integrate(prior, prior.pivots .* 0 .+ 1) ≈ 1

@test mean(prior) ≈ .4

cprior = condition(prior, low = .3)
@test bad.integrate(cprior, prior.pivots .* 0 .+ 1) ≈ 1
@test mean(cprior) > mean(prior)
