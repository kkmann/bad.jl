using Test

mle = MaximumLikelihoodEstimator()
rbe = RaoBlackwellEstimator()

x1_early_stop = 0:n1(design)
x1_early_stop = x1_early_stop[n2.(design, x1_early_stop) .== 0]

@test all( mle.(x1_early_stop, 0, design) .== rbe.(x1_early_stop, 0, design) )

@test bias.(0:.01:1, rbe, design) |> b -> abs.(b) |> maximum < 1e-12

XX = sample_space(design)
@test all(0 .<= rbe.(XX[:,1], XX[:,2], design) .<= 1)
