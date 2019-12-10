XX = sample_space(design)

@test 1 ≈ pmf.(design, XX[:,1], XX[:,2], pmcr) |> sum
@test 1 ≈ pmf.(design, XX[:,1], XX[:,2], prior) |> sum

X1 = 0:n1(design)
@test 1 ≈ pmf.(design, X1, pmcr) |> sum
@test 1 ≈ pmf.(design, X1, prior) |> sum

for x1 in X1
    X2 = sample_space(design, x1)
    @test 1 ≈ pmf_x2_given_x1.(design, x1, X2, pmcr) |> sum
    @test 1 ≈ pmf_x2_given_x1.(design, x1, X2, prior) |> sum
end
