using Test

pnull = .2
pmcr  = .3
prior = Beta(5, 7)
α, β  = .05, .2

samplesize = SampleSize(prior)

maximise(samplesize)

obj   = minimise(samplesize)

pwr   = subject_to(Power(prior, pmcr), β)
toer  = subject_to(TypeOneErrorRate(PointMass(pnull), pnull), α)

design = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(pnull, α),
        minimal_expected_power(prior, pmcr, 1 - β),
    ) |> optimise

@test obj(design) == design.score
@test pwr(design) >= 1 - pwr.β
@test toer(design) <= toer.α
