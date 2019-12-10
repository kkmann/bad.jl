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
@test !mlecompatible(design, p0, α)["compatible"]

cmle = CompatibleMLE(design)
@test compatible(EstimatorOrdering(cmle), design, p0, α)["compatible"]

XX    = sample_space(design)
resid = cmle.(XX[:,1], XX[:,2], design) .- mle.(XX[:,1], XX[:,2], design)
@test .001 < maximum(resid)



α, β = .05, .2
function designh0(p0, p1)
    mtoer = TypeOneErrorRate(PointMass(p0), p0)
    power = Power(PointMass(p1), p1)
    problem = Problem(
        minimise(SampleSize(PointMass(p0))),
        subject_to(mtoer, α),
        subject_to(power, β)
    )
    design = optimise(problem)
end
function designh1(p0, p1)
    mtoer = TypeOneErrorRate(PointMass(p0), p0)
    power = Power(PointMass(p1), p1)
    problem = Problem(
        minimise(SampleSize(PointMass(p1))),
        subject_to(mtoer, α),
        subject_to(power, β)
    )
    design = optimise(problem)
end

function test_design(design, p0)

    cmle  = CompatibleMLE(design)
    @test compatible(EstimatorOrdering(cmle), design, p0, α)["compatible"]
    mle   = MaximumLikelihoodEstimator()
    XX    = sample_space(design)
    resid = cmle.(XX[:,1], XX[:,2], design) .- mle.(XX[:,1], XX[:,2], design)
    if mlecompatible(design, p0, α)["compatible"]
        @test maximum(resid) < 5*1e-4 # numerical inaccuracy
    else
        @test maximum(resid) < 5*1e-2
    end
end

for p0 = 0.1:.1:.7
        p1 = p0 .+ .2
        test_design(designh0(p0, p1), p0)
        test_design(designh1(p0, p1), p0)
end
