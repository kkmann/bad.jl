using Test



function simon(r1::Int, n1::Int, r::Int, n::Int)
    nn2 = zeros(Int, n1 + 1)
    nn2[(r1 + 2):end] .= n - n1
    cc = repeat([Inf], n1 + 1)
    cc[(r1 + 2):end] .= r
    cc2 = cc .- collect(0:n1)
    cc2[cc2 .< 0]     .= -Inf
    cc2[cc2 .>= nn2]  .= Inf
    nn2[cc2 .== -Inf] .= 0
    Design(nn2, cc2)
end



α, β = 0.05, .2



p0, p1  = .05, .25
design1 = simon(0, 9, 2, 17)
design2 = Problem(
    minimise(SampleSize(PointMass(p0))),
    subject_to(TypeOneErrorRate(PointMass(p0), p0), α, [-.1, 1.1]),
    subject_to(Power(PointMass(p1), p1), β, [-.1, 1.1]),
    nmax                    = 50,
    n1values                = collect(5:10),
    n2maxcontinuereln1      = 6.,
    n2mincontinuereln1      = 0.0,
    n2mincontinueabs        = 10,
    type                    = :GroupSequential
)  |> pr -> optimise(pr, timelimit = 300)
pow     = Power(PointMass(p1), p1)
@test all( pow.([design1, design2]) .> 0.8)
err     = TypeOneErrorRate(PointMass(p0), p0)
@test all( err.([design1, design2]) .< 0.05)
ess     = SampleSize(PointMass(p0))
@test ess(design1) >= ess(design2)



p0, p1  = .1, .3
design1 = simon(1, 10, 5, 29)
design2 = Problem(
    minimise(SampleSize(PointMass(p0))),
    subject_to(TypeOneErrorRate(PointMass(p0), p0), α, [-.1, 1.1]),
    subject_to(Power(PointMass(p1), p1), β, [-.1, 1.1]),
    nmax                    = 50,
    n1values                = collect(5:15),
    n2maxcontinuereln1      = 7.,
    n2mincontinuereln1      = 0.0,
    n2mincontinueabs        = 7,
    curtail_stage_one_fct   = 0.0,
    type                    = :GroupSequential
) |> pr -> optimise(pr, timelimit = 300)
pow     = Power(PointMass(p1), p1)
@test all( pow.([design1, design2]) .> 0.8)
err     = TypeOneErrorRate(PointMass(p0), p0)
@test all( err.([design1, design2]) .< 0.05)
ess     = SampleSize(PointMass(p0))
@test ess(design1) >= ess(design2)



p0, p1  = .7, .9
design1 = simon(4, 6, 22, 27)
design2 = Problem(
    minimise(SampleSize(PointMass(p0))),
    subject_to(TypeOneErrorRate(PointMass(p0), p0), α, [-.1, 1.1]),
    subject_to(Power(PointMass(p1), p1), β, [-.1, 1.1]),
    nmax                    = 50,
    n1values                = collect(5:10),
    n2maxcontinuereln1      = 6.,
    n2mincontinuereln1      = 0.0,
    n2mincontinueabs        = 10,
    type                    = :GroupSequential
) |> pr -> optimise(pr, timelimit = 300)
pow     = Power(PointMass(p1), p1)
@test all( pow.([design1, design2]) .> 0.8)
err     = TypeOneErrorRate(PointMass(p0), p0)
@test all( err.([design1, design2]) .< 0.05)
ess     = SampleSize(PointMass(p0))
@test ess(design1) >= ess(design2)
