prior = Beta(mean = .35, sd = .1)
p0    = .2
α, β  = .05, .2

get_optimal_design(prior, p0, α, β, type) =
    Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .05, 1 - β),
        type = type
    )

ts  = get_optimal_design(prior, p0, α, β, :TwoStage) |> optimise
gs  = get_optimal_design(prior, p0, α, β, :GroupSequential) |> optimise
os  = get_optimal_design(prior, p0, α, β, :OneStage) |> optimise
