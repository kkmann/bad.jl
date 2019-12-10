using Test, bad; import Plots, Distributions

p = 0:0.1:1
for n in 0:15, pp in p
    x = 0:n
    @test all( Distributions.pdf.(Distributions.Binomial(n, pp), x) .≈ bad.dbinom.(x, n, pp) )
    @test all( Distributions.cdf.(Distributions.Binomial(n, pp), x) .≈ bad.pbinom.(x, n, pp) )
end

for a in .5:.5:10, b in .5:.5:10
    @test all( Distributions.pdf.(Distributions.Beta(a, b), p) .≈ bad.dbeta.(p, a, b) )
    @test all( Distributions.cdf.(Distributions.Beta(a, b), p) .≈ bad.pbeta.(p, a, b) )
end

@test true
