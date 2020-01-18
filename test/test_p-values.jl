using Test, bad


p0  = .2
α = .05
p1 = .4
β = .2
mcr = .3
pmp = .7

prior = (0.2*Beta(1, 1) + 0.8*Beta(mean = .35, sd = .1)) <= pmp

problem = Problem(
    minimise(
        SampleSize(prior)
    ),
    Power(prior  | p0) <= α,
    Power(prior >= mcr)  >= 1 - β
)
design = optimise(problem; verbosity = 0)



𝚾 = sample_space(design)

ordering = EstimatorOrdering(MaximumLikelihoodEstimator(); orientation = :superiority)
pval     = PValue(ordering, design, p0)
pvals    = pval.(𝚾[:,1], 𝚾[:,2])

@test all(
    pvals .≈ p_value.(𝚾[:,1], 𝚾[:,2], p0, ordering, design)
)

prob0 = pmf.(𝚾[:,2], n2.(design, 𝚾[:,1]), p0) .* pmf.(𝚾[:,1], n1(design), p0)

function f(p)
    inds = pvals .<= p
    min(1, max(0, sum(prob0[inds])))
end

p = 0:.01:1
@test all(f.(p) .<= p)
