function simon(p0, p1, α, β)
    Problem(
        minimise_expected_sample_size(PointMass(p0)),
        maximal_type_one_error_rate(p0, α, k = 5),
        minimal_expected_power(PointMass(p1), p1, 1 - β,
        conditional_threshold = -1, power_curtail = 2),
        type = :GroupSequential,
        max_rel_increase = 5,
        min_rel_increase = 1
    ) |>
    x -> optimise(x, timelimit = 900)
end

function simon(r1::Int, n1::Int, r::Int, n::Int)
    nn2 = zeros(Int, n1 + 1)
    nn2[(r1 + 2):end] .= n - n1
    cc = repeat([Inf], n1 + 1)
    cc[(r1 + 2):end] .= r
    cc2 = cc .- collect(0:n1)
    Design(nn2, cc2)
end

α, β = .05, .2

p0, p1  = .05, .25
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(0, 9, 2, 17)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)

p0, p1  = .1, .3
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(1, 10, 5, 29)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)
@test n1(design2)  == n1(design1)

p0, p1  = .2, .4
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(3, 13, 12, 43)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)
@test n1(design2)  == n1(design1)

p0, p1  = .3, .5
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(5, 15, 18, 46)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)
@test n1(design2)  == n1(design1)

p0, p1  = .4, .6
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(7, 16, 23, 46)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)
@test n1(design2)  == n1(design1)

p0, p1  = .5, .7
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(8, 15, 26, 43)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)
@test n1(design2)  == n1(design1)

p0, p1  = .6, .8
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(7, 11, 30, 43)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)

p0, p1  = .7, .9
ess     = bad.ExpectedSampleSize(PointMass(p0))
design1 = simon(4, 6, 22, 27)
design2 = simon(p0, p1, α, β)
@test ess(design2) <= ess(design1)
