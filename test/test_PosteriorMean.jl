using Test

𝚾 = sample_space(design)

pme1 = PosteriorMean(prior)
pme2 = PosteriorMeanPrecalculated(prior, design)

@test all( pme1.(𝚾[:,1], 𝚾[:,2], design) ≈ pme2.(𝚾[:,1], 𝚾[:,2], design) )
@time pme1.(𝚾[:,1], 𝚾[:,2], design)
@time pme2.(𝚾[:,1], 𝚾[:,2], design)

pme3 = PosteriorMean(JeffreysPrior(design))
pme4 = PosteriorMeanPrecalculated(JeffreysPrior(design), design)
@test all( pme3.(𝚾[:,1], 𝚾[:,2], design) ≈ pme4.(𝚾[:,1], 𝚾[:,2], design) )
@time pme3.(𝚾[:,1], 𝚾[:,2], design)
@time pme4.(𝚾[:,1], 𝚾[:,2], design)

import Plots
p   = 0:.01:1

Plots.plot(p, hcat(
        bias.(p, pme2, design),
        bias.(p, pme4, design)
))


Plots.plot(p, hcat(
        sqrt.(mean_squared_error.(p, pme2, design) ),
        sqrt.(mean_squared_error.(p, pme4, design) )
))
