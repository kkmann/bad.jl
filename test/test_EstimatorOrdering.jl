p0, α     = .2, .05
prior, β  = PointMass(0.4), .2
design    = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .05, 1 - β),
) |> optimise

mle      = MaximumLikelihoodEstimator()
ordering = EstimatorOrdering(mle)
space    = sample_space(design)

smaller_or_equal(space[1, 1], space[1, 2], space[2, 1], space[2, 2], ordering, design)

p_value(5, 2, p0, ordering, design; orientation = :superiority)

ci = ClopperPearsonInterval(ordering, design)

ci(15, 0)
