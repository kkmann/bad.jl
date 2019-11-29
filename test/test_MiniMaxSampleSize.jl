prior = PointMass(.4)
p0    = .2
α, β  = .05, .2

obj  = MiniMaxSampleSize()
toer = maximal_type_one_error_rate(p0, α)
pwr  = minimal_expected_power(prior, p0 + .2, 1 - β)

design = Problem(obj, toer, pwr) |> optimise

plot(design)
