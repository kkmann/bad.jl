prior = Beta(mean = .5, sd = .05)
p0    = .2
α, β  = .05, .2

get_optimal_design(prior, p0, α, β; group_sequential = false, one_stage = false) =
    DesignIPModel(prior, p0, α, β; group_sequential = group_sequential, one_stage = one_stage) +
    minimal_expected_power(prior, p0, 1 - β) +
    minimize_expected_sample_size() |>
    optimise

ts  = get_optimal_design(prior, p0, α, β)
gs  = get_optimal_design(prior, p0, α, β; group_sequential = true)
os  = get_optimal_design(prior, p0, α, β; one_stage = true)
os2 = DesignIPModel(prior, p0, α, β;
        nmax = 32, n1min = 32, n1max = 32,
        one_stage = true
    ) +
    minimal_expected_power(prior, p0, 1 - β) +
    minimize_expected_sample_size() |>
    optimise
