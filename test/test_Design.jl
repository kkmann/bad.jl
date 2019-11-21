prior  = Beta(mean = .4, sd = .1)
p0     = .2
α, β   = .05, .2
design = get_optimal_design(prior, p0, α, β; verbose = 3)
