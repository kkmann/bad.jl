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

XX    = sample_space(design)
resid = cmle.(XX[:,1], XX[:,2], design) .- mle.(XX[:,1], XX[:,2], design)
@test .001 < maximum(resid)



α, β = .05, .2
design_h0 = (p0, p1) -> Problem(
                minimise_expected_sample_size(PointMass(p0)),
                maximal_type_one_error_rate(p0, α),
                minimal_expected_power(PointMass(p1), p0, 1 - β),
        ) |> pr -> optimise(pr, mle_incompatible_is_error = false)
design_h1 = (p0, p1) -> Problem(
                minimise_expected_sample_size(PointMass(p1)),
                maximal_type_one_error_rate(p0, α),
                minimal_expected_power(PointMass(p1), p0, 1 - β),
        ) |> pr -> optimise(pr, mle_incompatible_is_error = false)

function test_design(design, p0)
    cmle  = CompatibleMLE(design)
    @test compatible(cmle, design, p0, α)[1]
    mle   = MaximumLikelihoodEstimator()
    XX    = sample_space(design)
    resid = cmle.(XX[:,1], XX[:,2], design) .- mle.(XX[:,1], XX[:,2], design)
    if compatible(mle, design, p0, α)[1]
        @test maximum(resid) < 1e-4 # numerical inaccuracy
    else
        @test maximum(resid) < 1e-2
    end
end

for p0 = 0.1:.1:.7
        p1 = p0 .+ .2
        test_design(design_h0(p0, p1), p0)
        test_design(design_h1(p0, p1), p0)
end
