prior = Beta(mean = .4, sd = .1)
p0    = .2
α, β  = .05, .2

obj  = minimise_expected_sample_size(prior)
toer = maximal_type_one_error_rate(p0, α)
pwr  = minimal_expected_power(prior, p0 + .1, 1 - β)

get_optimal_design(prior, p0, α, β, type) = Problem(obj, toer, pwr, type = type)

ts  = get_optimal_design(prior, p0, α, β, :TwoStage) |> optimise
gs  = get_optimal_design(prior, p0, α, β, :GroupSequential) |> optimise
os  = get_optimal_design(prior, p0, α, β, :OneStage) |> optimise

tbl = as_table(ts)

plot(ts)

@test obj(ts) < obj(gs) < obj(os)
for design in (ts, gs, os)
    @test obj(design) ≈ design.score
    @test pwr(design) .>= 0
    @test toer(design) .<= 0
    @test all(pwr.(design, (early_futility(design) + 1):n1(design)) .>= 0)
end


# prior = Beta(mean = .35, sd = .1)
# p0    = .2
# α, β  = .05, .2
#
# get_optimal_design(prior, p0, α, β, type) =
#     Problem(
#         minimise_expected_sample_size(prior),
#         maximal_type_one_error_rate(p0, α),
#         minimal_expected_power(prior, p0 + .05, 1 - β),
#         n1min = 25,
#         n1max = 150,
#         nmax  = 150
#     )
#
# ts  = get_optimal_design(prior, p0, α, β, :TwoStage) |> p _> optimise
