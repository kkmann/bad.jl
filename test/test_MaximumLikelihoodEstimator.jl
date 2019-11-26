p0, α     = .2, .05
prior, β  = PointMass(0.4), .1
ts        = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .05, 1 - β),
) |> optimise

mle = MaximumLikelihoodEstimator()
p   = 0:.01:1
bias.(p, mle, ts)
sqrt.(mean_squared_error.(p, mle, ts))



space = sample_space(ts)
probability.(space[:,1], space[:,2], ts, .3) |> sum
