using Test

jprior = JeffreysPrior(design)

p = .025:.01:.975
@test all( pdf.(p, jprior) .> 0 )

posterior = update(jprior, 7, 12)

@test mean(posterior) > mean(jprior)
