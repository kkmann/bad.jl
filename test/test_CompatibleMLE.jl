using Test

α, β     = .05, .1
p0       = .6
shan_n1  = 19
shan_n2  = zeros(19 + 1)
shan_n2[(13:17) .+ 1] .= [33, 31, 31, 31, 14]
shan_c   = vcat(
        [Inf for x1 in 0:12],
        [36, 35, 35, 35, 22],
        [-Inf for x1 in 18:shan_n1]
)
design = Design(shan_n2, shan_c .- (0:shan_n1))

mle  = MaximumLikelihoodEstimator()
@test !compatible(mle, design, p0, α)[1]

cmle = CompatibleMLE(design)
@test compatible(cmle, design, p0, α)[1]



p0     = .2
α, β   = .05, .2
prior  = PointMass(0.4)
design = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .05, 1 - β),
) |> optimise

mle  = MaximumLikelihoodEstimator()
@test compatible(mle, design, p0, α)[1]

cmle = CompatibleMLE(design)
@test compatible(cmle, design, p0, α)[1]

XX = sample_space(design)

@test abs.(
        cmle.(XX[:,1], XX[:,2], design) .-
        mle.(XX[:,1], XX[:,2], design)
) |> maximum < 1e-4
