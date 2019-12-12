using Test, bad; import Plots


pnull = .25
palt  = .4
p     = Beta(1, 1)
α, β  = .05, .2
problem = Problem(
    minimise(SampleSize(p | palt)),
    subject_to(TypeOneErrorRate(p | pnull), α),
    subject_to(Power(p | palt), β)
)
design = optimise(problem; verbosity = 0)

samplesize = SampleSize(p)

@test n1(design) == samplesize(design, 0, 0, .2)
@test n1(design) ≈ samplesize(design, 0, 0)
@test all( n1(design) .≈ samplesize.(design, early_stop_region(design), .2) )
@test all( n1(design) .≈ samplesize.(design,  early_stop_region(design)) )
@test all( n1(design) .< samplesize.(design,  continuation_region(design)) )



pow = problem.power.score

@test 0 ≈ pow(design, 0, 0, .2)
@test (0 .≈ pow.(design, futility_region(design), .4) ) |> all
@test (0 .≈ pow.(design, futility_region(design)) ) |> all
@test (1 .≈ pow.(design, efficacy_region(design), .4) ) |> all
@test (1 .≈ pow.(design, efficacy_region(design)) ) |> all
@test (0 .≈ pow.(design, continuation_region(design), 0)) |> all
@test (0.5 .< pow.(design, continuation_region(design)) .< .99) |> all
@test 0 ≈ pow(design, 6, .2)
@test 0 ≈ pow(design, .2)
@test abs(1 - β - pow(design)) < 1e-3



toer = problem.toer.score

@test 0 ≈ toer(design, 0, 0, .2)
@test (0 .≈ toer.(design, futility_region(design), .1) ) |> all
@test (0 .≈ toer.(design, futility_region(design)) ) |> all
@test (1 .≈ toer.(design, efficacy_region(design), .1) ) |> all
@test (1 .≈ toer.(design, efficacy_region(design)) ) |> all
@test (0 .≈ toer.(design, continuation_region(design), 0)) |> all
@test (0.0 .< toer.(design, continuation_region(design)) .< .99) |> all
@test 1 ≈ toer(design, 20, pnull)
@test toer(design, pnull) <= α
@test toer(design, pnull/2) < toer(design, pnull)

# check for monotone conditional error/power (for n1partial > 0)
for n1partial in 1:n1(design)
    println(n1partial)
    @test begin
        map(
            x1p -> power(design; x1partial = x1p, n1partial = n1partial),
            collect(0:n1partial)
        ) |>
        x -> diff(x) |>
        x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
    end
    @test begin
        map(
            x1p -> toer(design; x1partial = x1p, n1partial = n1partial),
            collect(0:n1partial)
        ) |>
        x -> diff(x) |>
        x -> minimum(x) >= -sqrt(eps()) # numerical inaccuracies
    end
end
