p0, α     = .2, .05
α, β      = .05, .2
prior     = condition(.2*Beta(1, 1) + .8*Beta(5, 7), high = .6)
design    = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .1, 1 - β),
) |> optimise

mle  = MaximumLikelihoodEstimator()
pme1 = PosteriorMean(prior)
pme2 = PosteriorMean(JeffreysPrior(design))

using Plots
p   = 0:.01:1

Plots.plot(p, hcat(
        bias.(p, mle, design),
        bias.(p, pme1, design),
        bias.(p, pme2, design),
))


Plots.plot(p, hcat(
        sqrt.(mean_squared_error.(p, mle, design) ),
        sqrt.(mean_squared_error.(p, pme1, design) ),
        sqrt.(mean_squared_error.(p, pme2, design) ),
))
