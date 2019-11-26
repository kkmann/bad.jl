prior = Beta(mean = .4, sd = .1)
p0    = .2
α, β  = .05, .2

problem = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .1, 1 - β)
    )

design1 = optimise(problem)

design1.model

xx1, nn1 = 3, 7

design2 = adapt(design1, xx1, nn1)
