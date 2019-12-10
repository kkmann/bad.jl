using Test, bad; import Plots


prior1 = Beta(mean = .4, sd = .1)
@test expectation(p -> 1, prior1) ≈ 1
@test (prior1.a ≈ 9.2) & (prior1.b ≈ 13.8)
@test mean(prior1) ≈ .4
prior_sd = expectation(p -> (p - mean(prior1))^2, prior1) |> sqrt
@test prior_sd ≈ .1

prior2 = condition(prior1; low = .3)
@test minimum(prior2.pivots) > .3
@test expectation(p -> 1, prior2) ≈ 1
@test (prior2.a ≈ 9.2) & (prior2.b ≈ 13.8)
@test (prior2.low ≈ .3) & (prior2.high ≈ 1.0)
@test mean(prior2) > mean(prior1)

prior3 = condition(prior1; high = .5)
@test maximum(prior3.pivots) < .5
@test expectation(p -> 1, prior3) ≈ 1
@test (prior3.a ≈ 9.2) & (prior3.b ≈ 13.8)
@test (prior3.low ≈ .0) & (prior3.high ≈ .5)
@test mean(prior3) < mean(prior1)

prior4 = Beta(2, 3)
@test expectation(p -> 1, prior4) ≈ 1
@test (prior4.a ≈ 2) & (prior4.b ≈ 3)
@test (prior4.low ≈ .0) & (prior4.high ≈ 1.)
@test mean(prior4) ≈ 2/5

post1 = update(prior4, 1, 1)
@test expectation(p -> 1, post1) ≈ 1
@test (post1.a ≈ 3) & (post1.b ≈ 3)
@test (post1.low ≈ .0) & (post1.high ≈ 1.)

post2 = update(prior3, 0, 1)
@test expectation(p -> 1, post2) ≈ 1
@test (post2.a ≈ 9.2) & (post2.b ≈ 14.8)
@test (post2.low ≈ .0) & (post2.high ≈ .5)



# mixtures

prior1  = condition(Beta(mean = .4, sd = .1), low = .3)
prior2  = Beta(1, 1)

mprior1 = .2prior1 + .3prior2
@test !is_proper(mprior1)
p = .01:.01:.99
@test all( .2*pdf.(p, prior1) + .3*pdf.(p, prior2) .== pdf.(p, mprior1) )

mprior2    = .9prior1 + .1prior2
@test mean(prior1) < mean(mprior2) < mean(prior2)
prior1_sd  = expectation(p -> (p - mean(prior1))^2, prior1) |> sqrt
mprior2_sd = expectation(p -> (p - mean(mprior2))^2, mprior2) |> sqrt
@test prior1_sd < mprior2_sd

mprior3 = .5prior1 + .5prior2
@test mean(mprior2) < mean(mprior3) < mean(prior2)
prior2_sd  = expectation(p -> (p - mean(prior2))^2, prior2) |> sqrt
mprior3_sd = expectation(p -> (p - mean(mprior3))^2, mprior3) |> sqrt
@test prior2_sd > mprior3_sd

mpost = update(mprior2, 1, 1)
@test (mpost.priors[1].a ≈ 10.2) & (mpost.priors[1].b ≈ 13.8)
@test (mpost.priors[2].a ≈ 2) & (mpost.priors[2].b ≈ 1)
@test mean(mpost) > mean(mprior2)

mpost2 = update(mprior2, 2, 8)
@test mean(mpost2) < mean(mprior2)
mpost2_sd = expectation(p -> (p - mean(mpost2))^2, mpost2) |> sqrt
@test mpost2_sd < mprior2_sd
