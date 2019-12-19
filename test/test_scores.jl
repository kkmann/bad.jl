using Test, bad; import Plots, Distributions


pnull = .25
palt  = .4
p     = Beta(5, 7)
α, β  = .05, .2
problem = Problem(
        minimise( SampleSize(p) ),
        Power(p | pnull) <= α,
        Power(p | palt)  >= 1 - β,
        nmax  = 100,
        n1max = 40
    )
design = optimise(problem; verbosity = 0)



samplesize = SampleSize(p)

@test n1(design) == samplesize(design, 0, 0, .2)
@test n1(design) ≈ samplesize(design, 0, 0)
@test all( n1(design) .≈ samplesize.(design, early_stop_region(design), .2) )
@test all( n1(design) .≈ samplesize.(design,  early_stop_region(design)) )
@test all( n1(design) .< samplesize.(design,  continuation_region(design)) )
@test expectation( p -> sum(
        Distributions.pdf.(Distributions.Binomial(n1(design), p), 0:n1(design)) .*
        n.(design, 0:n1(design))
    ),
    p
) ≈ samplesize(design)


pow = problem.power.score

@test 0 ≈ pow(design, 0, 0, .2)
@test (0 .≈ pow.(design, futility_region(design), pow.bounds[1]) ) |> all
@test (0 .≈ pow.(design, futility_region(design)) ) |> all
@test (1 .≈ pow.(design, efficacy_region(design), pow.bounds[1]) ) |> all
@test (1 .≈ pow.(design, efficacy_region(design)) ) |> all
@test (0 .≈ pow.(design, continuation_region(design), 0.) ) |> all
@test (0.5 .< pow.(design, continuation_region(design)) .< .99) |> all
@test 0 ≈ pow(design, 6, .2)
@test 0 ≈ pow(design, .2)
@test abs(1 - β - pow(design)) < 1e-3



toer = problem.toer.score

@test 0 ≈ toer(design, 0, 0, .2)
@test (0 .≈ toer.(design, futility_region(design), toer.bounds[1]) ) |> all
@test (0 .≈ toer.(design, futility_region(design)) ) |> all
@test (1 .≈ toer.(design, efficacy_region(design), toer.bounds[1]) ) |> all
@test (1 .≈ toer.(design, efficacy_region(design)) ) |> all
@test (0 .≈ toer.(design, continuation_region(design), 0)) |> all
@test (0.0 .< toer.(design, continuation_region(design)) .< .99) |> all
@test 1 ≈ toer(design, 20, pnull)
@test toer(design, pnull) <= α
@test toer(design, pnull/2) < toer(design, pnull)

# check for monotone conditional error/power (for n1partial > 0)
for n1partial in 1:n1(design)
    @test begin
        map(
            x1p -> pow(design; partial_stage_one = (x1p, n1partial) ),
            0:n1partial
        ) |>
        x -> diff(x) |>
        x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
    end
    @test begin
        map(
            x1p -> toer(design; partial_stage_one = (x1p, n1partial) ),
            0:n1partial
        ) |>
        x -> diff(x) |>
        x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
    end
end

# check for monotone conditional error/power in stage two
for x1 in 0:n1(design)
    for n2partial in 1:n2(design, x1)
        @test begin
            map(
                x2p -> pow(design, x1; partial_stage_two = (x2p, n2partial) ),
                0:n2partial
            ) |>
            x -> diff(x) |>
            x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
        end
        @test begin
            map(
                x2p -> toer(design, x1; partial_stage_two = (x2p, n2partial) ),
                0:n2partial
            ) |>
            x -> diff(x) |>
            x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
        end
    end
end

@test true
