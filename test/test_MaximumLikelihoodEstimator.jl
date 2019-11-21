prior = PointMass(0.4)
p0    = .2
α, β  = .05, .2

ts = DesignIPModel(prior, p0, α, β) +
    minimal_expected_power(.8) +
    minimize_expected_sample_size() |>
    optimise

estimator = MaximumLikelihoodEstimator
