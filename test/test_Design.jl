sample_size    = [10, 10, 25, 23, 22, 10]
critical_value = [Futility(), Futility(), 11, 10, 9, Efficacy()]
stage_one_ss   = Int(5)
x1             = collect(0:stage_one_ss)

d = Design(sample_size, critical_value)

@test n1(d) == 5
@test n.(d, x1) == d.n
@test c.(d, x1) == d.c

@test_throws ErrorException n(d, -1)
@test_throws ErrorException n(d, stage_one_ss + 1)
@test_throws ErrorException c(d, -1)
@test_throws ErrorException c(d, stage_one_ss + 1)
