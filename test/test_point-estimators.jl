using Test, bad



pnull = .25
pmcr  = .3
p     = Beta(5, 7)
α, β  = .05, .2
problem = Problem(
    minimise(
        SampleSize(p | (pnull + .2) )
    ),
    Power(p  | pnull) <= α,
    Power(p >= pmcr)  >= 1 - β
)
design = optimise(problem; verbosity = 0)



𝚾 = sample_space(design)

pme1 = PosteriorMean(p)
pme2 = PosteriorMeanPrecalculated(p, design)

@test all( pme1.(𝚾[:,1], 𝚾[:,2], design) ≈ pme2.(𝚾[:,1], 𝚾[:,2], design) )

pme3 = PosteriorMean(JeffreysPrior(design))
pme4 = PosteriorMeanPrecalculated(JeffreysPrior(design), design)
@test all( pme3.(𝚾[:,1], 𝚾[:,2], design) ≈ pme4.(𝚾[:,1], 𝚾[:,2], design) )




mle = MaximumLikelihoodEstimator()
rbe = RaoBlackwellEstimator()

x1_early_stop = early_stop_region(design)
pp = 0:.01:1
@test all( mle.(x1_early_stop, 0, design) .== rbe.(x1_early_stop, 0, design) )
# design is injective on continuation region, rbe = x1/n1?
# allow numerical inaccuracy for rbe!
@test all( rbe.(𝚾[:,1], 𝚾[:,2], design) .≈ 𝚾[:,1]./n1(design) )
@test maximum( abs.( bias.(pp, rbe, design) ) ) ≈ 0 atol = 1e-12
@test all(0 .<= rbe.(𝚾[:,1], 𝚾[:,2], design) .<= 1)




# build shan design
α, β     = .05, .1
pnull       = .6
shan_n1  = 19
shan_n2  = zeros(19 + 1)
shan_n2[(13:17) .+ 1] .= [33, 31, 31, 31, 14]
shan_c   = vcat(
        [Inf for x1 in 0:12],
        [36, 35, 35, 35, 22],
        [-Inf for x1 in 18:shan_n1]
)
design = Design(shan_n2, shan_c .- (0:shan_n1))

@test !mlecompatible(design, pnull, α)["compatible"]

cmle  = CompatibleMLE(design)
@test compatible(EstimatorOrdering(cmle), design, pnull, α)["compatible"]

𝚾     = sample_space(design)
resid = cmle.(𝚾[:,1], 𝚾[:,2], design) .- mle.(𝚾[:,1], 𝚾[:,2], design)
# make sure there is a (minute) difference
@test .001 < maximum(resid)



α, β = .05, .2

get_design(p_sample_size, pnull, palt) = Problem(
        minimise( SampleSize(p | p_sample_size) ),
        Power(p  | pnull) <= α,
        Power(p >= palt)  >= 1 - β;
        maxmultipleonestage = 2.5
    ) |> problem -> optimise(problem; verbosity = 0)



for (pnull, palt) in [(pnull, pnull + .2) for pnull in 0.1:.1:.7]

    for psamplesize in (pnull, palt)
        design = get_design(psamplesize, pnull, palt)
        cmle   = CompatibleMLE(design; lambda = .25)
        if !compatible(EstimatorOrdering(cmle), design, pnull, α)["compatible"]
            global tmp = design, pnull, palt, psamplesize
        end
        @test compatible(EstimatorOrdering(cmle), design, pnull, α)["compatible"]
        𝚾     = sample_space(design)
        resid = cmle.(𝚾[:,1], 𝚾[:,2], design) .- mle.(𝚾[:,1], 𝚾[:,2], design)
        if mlecompatible(design, pnull, α)["compatible"]
            @test maximum(resid) ≈ 0 atol = 1e-2 # numerical inaccuracy
        else
            @test maximum(resid) < 5*1e-1 # shouldnt be tooo different
        end
    end
end

@test true
