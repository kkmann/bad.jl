using Test, bad



pnull = .2
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

ordering = EstimatorOrdering(MaximumLikelihoodEstimator(); orientation = :superiority)
pval     = PValue(ordering, design, pnull)

all(more_extreme.(
    pval.ordered_sample_space[2:end,1],
    pval.ordered_sample_space[2:end,2],
    pval.ordered_sample_space[1:(end-1),1],
    pval.ordered_sample_space[1:(end-1),2],
    ordering, design))

pval.(𝚾[:,1], 𝚾[:,2])

@test all(
    pval.(𝚾[:,1], 𝚾[:,2]) .≈ p_value.(𝚾[:,1], 𝚾[:,2], pnull, ordering, design)
)
