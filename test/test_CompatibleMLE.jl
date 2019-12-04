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
        ) |> optimise
design_h1 = (p0, p1) -> Problem(
                minimise_expected_sample_size(PointMass(p1)),
                maximal_type_one_error_rate(p0, α),
                minimal_expected_power(PointMass(p1), p0, 1 - β),
        ) |> optimise

function test_design(design, p0)
        mle  = MaximumLikelihoodEstimator()
        @test compatible(mle, design, p0, α)[1]
        cmle = CompatibleMLE(design, λ = 1.)
        @test compatible(cmle, design, p0, α)[1]
        XX    = sample_space(design)
        resid = cmle.(XX[:,1], XX[:,2], design) .- mle.(XX[:,1], XX[:,2], design)
        println(p0)
        @test maximum(resid) < 1e-4
end

for p0 = 0.1:.1:.7
        p1 = p0 .+ .2
        test_design(design_h0(p0, p1), p0)
        test_design(design_h1(p0, p1), p0)
end

# p = PValue(EstimatorOrdering(cmle), design, p0)

# XX_ordered = p.ordered_sample_space
# @test all( bad.more_extreme.(
#                 XX_ordered[2:end,1], XX_ordered[2:end,2],
#                 XX_ordered[1:(end - 1),1], XX_ordered[1:(end - 1),2],
#                 cmle,
#                 design
#         ) )
#
# @test 1e-9 > maximum(abs.(p.(XX[:,1], XX[:,2]) - p_value.(XX[:,1], XX[:,2], p0, cmle, design) ))
#
# @time p.(XX[:,1], XX[:,2])
# @time p_value.(XX[:,1], XX[:,2], p0, EstimatorOrdering(cmle), design)
