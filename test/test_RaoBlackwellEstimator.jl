using Test

p0     = .2
α, β   = .05, .2
prior  = PointMass(0.4)
design = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .05, 1 - β),
) |> optimise

mle = MaximumLikelihoodEstimator()
rbe = RaoBlackwellEstimator()

x1_early_stop = 0:n1(design)
x1_early_stop = x1_early_stop[n2.(design, x1_early_stop) .== 0]

@test all( mle.(x1_early_stop, 0, design) .== rbe.(x1_early_stop, 0, design) )

XX = sample_space(design)

@test all(0 .<= rbe.(XX[:,1], XX[:,2], design) .<= 1)

@test bias.(0:.01:1, rbe, design) |> b -> abs.(b) |> maximum < 1e-12
