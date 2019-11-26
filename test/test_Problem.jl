prior = Beta(mean = .4, sd = .1)
p0    = .2
α, β  = .05, .2

problem = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .1, 1 - β)
    )

design1 = optimise(problem)

update!(problem, update(prior, 2, 3))

design2 = optimise(problem)
