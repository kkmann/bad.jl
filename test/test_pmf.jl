using Test, bad

import Distributions
import Plots



p = 0:0.1:1
for n in 0:15, x = 0:n, pp in p
    @test Distributions.pdf.(Distributions.Binomial(n, pp), x) ≈ pmf(x, n, pp)
    @test Distributions.cdf.(Distributions.Binomial(n, pp), x) ≈ cdf(x, n, pp)
end



p = Beta(5, 7)

# check that marginal distribution (of x1) sums to 1
for xpartial in 0:10, npartial in xpartial:10
    pmff = pmf.(0:10, 10, p; xpartial = xpartial, npartial = npartial)
    @test sum(pmff) ≈ 1.0
    @test cumsum(pmff) ≈ cdf.(0:10, 10, p; xpartial = xpartial, npartial = npartial)
end



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

design = simon(4, 18, 10, 33)
XX     = sample_space(design)
x1, x2 = XX[:,1], XX[:,2]



# check that joint distribution sums to for all possible partial x1 observations
for x1partial in 0:n1(design)
    @test 1.0 ≈ pmf.(x1, x2, design, p; x1partial = x1partial, n1partial = n1(design)) |> sum
end

# check that marginal distribution (of x1) sums to 1
for x1partial in 0:n1(design), n1partial in x1partial:n1(design)
    pmff = pmf.(0:n1(design), n1(design), p; xpartial = x1partial, npartial = n1partial)
    @test sum(pmff) ≈ 1.0
    @test cumsum(pmff) ≈ cdf.(0:n1(design), n1(design), p; xpartial = x1partial, npartial = n1partial)
end

# check that conditional distribution (of x2 given x2) sums to 1
for x1 in 0:n1(design)
    pmff = pmf.(0:n2(design, x1), n2(design, x1), p)
    @test sum(pmff) ≈ 1.0
    @test all( cumsum(pmff) .≈ cdf.(0:n2(design, x1), n2(design, x1), p))
end

@test true
