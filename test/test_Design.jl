prior  = Beta(mean = .4, sd = .1)
p0     = .2
α, β   = .05, .2
ts = get_optimal_design(prior, p0, α, β; verbose = 3)
gs = get_optimal_design(prior, p0, α, β; verbose = 3, group_sequential = true)
os = get_optimal_design(prior, p0, α, β; verbose = 3, one_stage = true)
os = get_optimal_design(
    prior, p0, α, β; nmax = 52, n1min=52,
    n1max =52, max_rel_increase = Inf, min_conditional_power = 0.,
    min_rel_increase= 0, verbose = 3, single_stage = true
)
