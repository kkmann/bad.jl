p0, α     = .2, .05
prior, β  = PointMass(0.4), .1
ts        = DesignIPModel(prior, p0, α, β) +
        minimal_expected_power(prior, p0, 1 - β) +
        minimize_expected_sample_size() |>
        optimise

mle = MaximumLikelihoodEstimator()
p   = 0:.01:1
bias.(p, mle, ts)
sqrt.(mean_squared_error.(p, mle, ts))



space = sample_space(ts)
probability.(space[:,1], space[:,2], ts, .3) |> sum
