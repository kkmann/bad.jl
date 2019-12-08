Pkg.activate("."); using Revise; using Test
using bad

p0, α     = .2, .05
prior, β  = Beta(5, 7), .2
design    = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, p0 + .1, 1 - β),
    ) |> optimise

XX = sample_space(design)

X2 = sample_space(design, 3)

function pdf_k_less_n1(x1, x2, k, x_leq_k, design, prior)
    x   = x1 + x2
    return expectation(
        p -> binomial(n2(design, x1), x2) *
            binomial(n1(design) - k, x1 - x_leq_k) *
            p^(x - x_leq_k) *
            (1 - p)^(n(design, x1) - k - x + x_leq_k),
        prior
    )
end
function pdf_k_geq_n1(x1, x2, k, x_n1_to_k, design, prior)
    dn2 = n2(design, x1) - (k - n1(design))
    dn2 < 0 ? (return 0) : nothing
    dx2 = x2 - x_n1_to_k
    dx2 < 0 ? (return 0) : nothing
    res = expectation(
        p -> binomial(dn2, dx2) * p^dx2 * (1 - p)^(dn2 - dx2) *
            binomial(n1(design), x1) * p^x1 * (1 - p)^(n1(design) - x1),
        prior
    )
    println([x1, x2, res])
    res / sum(bad.dbinom.((early_futility(design) + 1):(early_efficacy(design) - 1), n1(design), p))
end

function pdf_k_greater_n1(x1, x2, τ, x_n1_to_τ, design, prior)
    τ <= n1(design) ? error("hi") : nothing
    τ  > n(design, x1) ? (return 0.0) : nothing
    n_geq_tau_preimage = (0:n1(design))[n.(design, 0:n1(design)) .>= τ]
    return expectation(
        p ->
            binomial(
                    n2(design, x1) - (τ - n1(design)),
                    x2 - x_n1_to_τ
                ) *
                p^(x2 - x_n1_to_τ) *
                (1 - p)^(n2(design, x1) - (τ - n1(design)) - (x2 - x_n1_to_τ)) *
            binomial(
                    n1(design),
                    x1
                ) *
                p^x1 *
                (1 - p)^(n1(design) - x1) /
            sum(
                binomial.(n1(design), n_geq_tau_preimage) .*
                    p.^n_geq_tau_preimage .*
                    (1 - p).^(n1(design) .- n_geq_tau_preimage)
            ),
        prior
    )
end

x1, x2 = , 13
k, x_n1_to_k= n1(design) + 20, 1
dn2 = n2(design, x1) - (k - n1(design))
dx2 = x2 - x_n1_to_k
p = .4
binomial(dn2, dx2) * p^(dx2) * (1 - p)^(dn2 - dx2) *
binomial(n1(design), x1) * p^x1 * (1 - p)^(n1(design) - x1)

x1, x2 = XX[5,:]
n1(design)
n2(design, x1)
pdf_k_less_n1.(XX[:,1], XX[:,2], n1(design), 2, design, PointMass(.4)) |> sum

XX_x1 = sample_space(design, 3)

pdf_k_greater_n1.(XX[:,1], XX[:,2], n1(design) + 1, 0, design, PointMass(.4)) |>
    sum

pdf_k_geq_n1.(XX[:,1], XX[:,2], n1(design) + 1, 0, design, PointMass(.4)) |>
    sum

design

bad.dbinom.(3:5, n1(design), .4) |> sum
