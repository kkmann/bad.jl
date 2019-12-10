using Test



samplesize = SampleSize(prior)

@test n1(design) == samplesize(design, 0, 0, .2)
@test n1(design) ≈ samplesize(design, 0, 0)
@test all( n1(design) .≈ samplesize.(design, early_stop_region(design), .2) )
@test all( n1(design) .≈ samplesize.(design,  early_stop_region(design)) )
@test all( n1(design) .< samplesize.(design,  continuation_region(design)) )
@test abs(28.4 - samplesize(design, .2)) < 1e-1



pow = Power(prior, pmcr)

@test 0 ≈ pow(design, 0, 0, .2)
@test (0 .≈ pow.(design, futility_region(design), .4) ) |> all
@test (0 .≈ pow.(design, futility_region(design)) ) |> all
@test (1 .≈ pow.(design, efficacy_region(design), .4) ) |> all
@test (1 .≈ pow.(design, efficacy_region(design)) ) |> all
@test (0 .≈ pow.(design, continuation_region(design), 0)) |> all
@test (0.5 .< pow.(design, continuation_region(design)) .< .99) |> all
@test 0 ≈ pow(design, 6, .2)
@test 0 ≈ pow(design, .2)
@test α <= pow(design, pmcr) <= 1 - β
@test abs(1 - β - pow(design)) < 1e-3



toer = TypeOneErrorRate(prior, pnull)

@test 0 ≈ toer(design, 0, 0, .2)
@test (0 .≈ toer.(design, futility_region(design), .1) ) |> all
@test (0 .≈ toer.(design, futility_region(design)) ) |> all
@test (1 .≈ toer.(design, efficacy_region(design), .1) ) |> all
@test (1 .≈ toer.(design, efficacy_region(design)) ) |> all
@test (0 .≈ toer.(design, continuation_region(design), 0)) |> all
@test (0.0 .< toer.(design, continuation_region(design)) .< .99) |> all
@test 1 ≈ toer(design, 9, pnull)
@test toer(design, pnull) <= α
@test toer(design, pnull/2) < toer(design, pnull)
