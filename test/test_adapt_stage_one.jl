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

ts.model

function adapt(problem, xx1, nn1, old_design)
    m, ind, n1_selected = bad.get_IP_model(problem)
    add!(m, ind, problem.toer, problem, xx1, nn1, old_design)
    add!(m, ind, problem.power, problem, xx1, nn1, old_design)
end

m, ind, n1_selected = bad.get_IP_model(ts.model)
m
bad.add!(m, ind, ts.model.toer, ts.model, 10, 20, ts)
m
bad.add!(m, ind, ts.model.power, ts.model, 10, 20, ts)
m
