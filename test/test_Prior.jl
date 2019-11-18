import bad.dbeta, bad.pbeta

priors = Prior(sin, low = 0., high = .5),
    condition(Prior(p -> 1), low = .5),
    BetaPrior(2., 2., low = .25, high = .75)

ϵ = sqrt(eps())

for pr in priors
    mid = pr.low/2 + pr.high/2
    @test all(isnan.((pr.pdf(pr.low - ϵ), pr.pdf(pr.high + ϵ))))
    @test all((pr.pdf(pr.low), pr.pdf(mid), pr.pdf(pr.high)) .>= 0)
    @test pr.cdf(pr.low)  ≈ 0
    @test pr.cdf(pr.high) ≈ 1
    @test 0 <= pr.cdf(mid) <= 1
end

prior        = BetaPrior(.5, .5)
p            = collect(0:.1:1)
p[1], p[end] = p[1] + sqrt(eps()), p[end] - sqrt(eps())
# small deviation due to excluding boundaries
@test all(0 .<= (prior.pdf.(p) .- dbeta.(p, .5, .5)) ./ dbeta.(p, .5, .5) .< .001)
@test all(0 .<= abs.(prior.cdf.(p[2:end]) .- pbeta.(p[2:end], .5, .5)) ./ pbeta.(p[2:end], .5, .5) .< .001)

prior  = Prior(p -> 1.)
cprior = condition(prior, high = .5)
@test cprior.cdf(.5) ≈ 1
cprior = condition(prior, low = .5)
@test cprior.cdf(.5) ≈ 0
