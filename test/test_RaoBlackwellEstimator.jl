p0, α    = .2, .05
prior, β = PointMass(0.4), .1
design   = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .05, 1 - β),
) |> optimise

mle = MaximumLikelihoodEstimator()
rbe = RaoBlackwellEstimator()

x1_early_stop = 0:n1(design)
x1_early_stop = x1_early_stop[n2.(design, x1_early_stop) .== 0]

@test all( mle.(x1_early_stop, 0, design) .== rbe.(x1_early_stop, 0, design) )
