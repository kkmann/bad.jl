using Test, bad; import Plots


pnull = .2
pmcr  = 0.3
palt  = .4
p     = Beta(1, 1)
α, β  = .05, .2
problem = Problem(
    minimise(SampleSize(p | palt)),
    subject_to(TypeOneErrorRate(p | pnull), α),
    subject_to(Power(p | palt), β)
)
design = optimise(problem; verbosity = 0)

XX = sample_space(design)

@test 1 ≈ pmf.(design, XX[:,1], XX[:,2], pmcr) |> sum
@test 1 ≈ pmf.(design, XX[:,1], XX[:,2], p) |> sum

X1 = 0:n1(design)
@test 1 ≈ pmf.(design, X1, pmcr) |> sum
@test 1 ≈ pmf.(design, X1, p) |> sum

for x1 in X1
    X2 = sample_space(design, x1)
    @test 1 ≈ pmf_x2_given_x1.(design, x1, X2, pmcr) |> sum
    @test 1 ≈ pmf_x2_given_x1.(design, x1, X2, p) |> sum
end

@test true
