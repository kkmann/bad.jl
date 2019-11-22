p0, α     = .2, .05
prior, β  = PointMass(0.4), .2
design    = DesignIPModel(prior, p0, α, β) +
        minimal_expected_power(1 - β) +
        minimize_expected_sample_size() |>
        optimise

mle     = MaximumLikelihoodEstimator()
smaller = EstimatorOrdering(mle)
space   = sample_space(design)

smaller(space[1, 1], space[1, 2], space[2, 1], space[2, 2], design)

p_value(5, 2, p0, smaller, design)
