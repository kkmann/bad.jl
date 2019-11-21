prior1  = condition(Beta(mean = .4, sd = .1), low = .3)
prior2  = Beta(1, 1)

mprior1 = .2prior1 + .3prior2
@test !is_proper(mprior3)

mprior2    = .9prior1 + .1prior2
@test mean(prior1) < mean(mprior2) < mean(prior2)
prior1_sd  = expected_value(p -> (p - mean(prior1))^2, prior1) |> sqrt
mprior2_sd = expected_value(p -> (p - mean(mprior2))^2, mprior2) |> sqrt
@test prior1_sd < mprior2_sd

mprior3 = .5prior1 + .5prior2
@test mean(mprior2) < mean(mprior3) < mean(prior2)
prior2_sd  = expected_value(p -> (p - mean(prior2))^2, prior2) |> sqrt
mprior3_sd = expected_value(p -> (p - mean(mprior3))^2, mprior3) |> sqrt
@test prior2_sd > mprior3_sd

mpost = update(mprior2, 1, 1)
@test (mpost.priors[1].a ≈ 10.2) & (mpost.priors[1].b ≈ 13.8)
@test (mpost.priors[2].a ≈ 2) & (mpost.priors[2].b ≈ 1)
@test mean(mpost) > mean(mprior2)

mpost2 = update(mprior2, 2, 8)
@test mean(mpost2) < mean(mprior2)
mpost2_sd = expected_value(p -> (p - mean(mpost2))^2, mpost2) |> sqrt
@test mpost2_sd < mprior2_sd
