prior = PointMass(.4)
p0    = .2
α, β  = .05, .2

obj  = MiniMaxSampleSize()
toer = maximal_type_one_error_rate(p0, α; k = 5)
pwr  = minimal_expected_power(prior, p0 + .2, 1 - β, power_curtail = .0)

design = Problem(obj, toer, pwr, n1max = 40, nmax = 40,
    min_rel_increase = .9, max_rel_increase = 5, min_abs_increase = 0) |> optimise

design.model

plot(design)

n(design)

design = Problem(obj, toer, pwr, n1max = 40, nmax = 40, type = :OneStage) |> optimise

design
