prior  = Beta(mean = .4, sd = .1)
p0     = .2
α, β   = .05, .2
design = get_optimal_design(prior, p0, α, β; verbose = 3)
design = get_optimal_design(prior, p0, α, β; verbose = 3, group_sequential = true)

design = get_optimal_design(
    prior, p0, α, β; nmax = 53, n1min=52,
    n1max =52, max_rel_increase = Inf, min_conditional_power = 0.,
    min_rel_increase= 0, verbose = 3, single_stage = true
)
