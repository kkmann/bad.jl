sample_size    = [10, 10, 25, 23, 22, 10]
critical_value = [Futility(), Futility(), 11, 10, 9, Efficacy()]
stage_one_ss   = Int(5)
x1             = collect(0:stage_one_ss)
d              = Design(sample_size, critical_value)

@test n1(d)     == 5
@test n.(d, x1) == d.n
@test c.(d, x1) == d.c

@test all(.!valid.(d, (-1, stage_one_ss + 1)))
@test all(valid.(d, x1))

@test probability(3, 10, d, .4) ≈ 7.20513705227343e-10

@test power(3, d, .4) ≈ 0.4365591524965087
@test power(d, .4) ≈ 0.2519087750133824

prior = Prior(mean = .4, sd = .1)
