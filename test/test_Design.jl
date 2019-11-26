prior = Beta(mean = .35, sd = .1)
p0    = .2
α, β  = .05, .2

get_optimal_design(prior, p0, α, β, type) =
    Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .05, 1 - β; conditional_threshold = .0, power_curtail = .999),
        type = type
    ) |>
    optimise

ts  = get_optimal_design(prior, p0, α, β, :TwoStage)
gs  = get_optimal_design(prior, p0, α, β, :GroupSequential)
os  = get_optimal_design(prior, p0, α, β, :OneStage)
os2 = DesignIPModel(prior, p0, α, β;
        nmax = 32, n1min = 32, n1max = 32,
        one_stage = true
    ) +
    minimal_expected_power(prior, p0, 1 - β) +
    minimize_expected_sample_size() |>
    optimise
