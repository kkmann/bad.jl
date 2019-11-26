prior = Beta(mean = .35, sd = .05)
p0    = .2
α, β  = .05, .2

mean(condition(prior, low = p0))

problem = Problem(
    minimise_expected_sample_size(prior),
    maximal_type_one_error_rate(p0, α),
    minimal_expected_power(prior, p0, 1 - β)
)

design = optimise(problem)

prior = Beta(mean = .35, sd = .05)
p0    = .2
α, β  = .05, .2

problem = Problem(
    minimise_expected_sample_size(prior),
    maximal_type_one_error_rate(p0, α),
    minimal_expected_power(prior, p0, 1 - β)
)

bad.n1vals(problem)

design = optimise(problem)

prior = Beta(mean = .35, sd = .1)
p0    = .2
α, β  = .05, .2

problem = Problem(
    minimise_expected_sample_size(prior),
    maximal_type_one_error_rate(p0, α),
    minimal_expected_power(prior, .3, 1 - β)
)

bad.n1vals(problem)

design = optimise(problem)
