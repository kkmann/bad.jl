using Test

p0, p1 = .2, .4
α, β   = .05, .2

# Simons minimax design (Optimal Two-Stage Designs for Phase II Clinical Trials, Richard Simon)
# r1/n1 = 4/18; r/n = 10/33
function simon(r1::Int, n1::Int, r::Int, n::Int)
    nn2 = zeros(Int, n1 + 1)
    nn2[(r1 + 2):end] .= n - n1
    cc = repeat([Inf], n1 + 1)
    cc[(r1 + 2):end] .= r
    cc2 = cc .- collect(0:n1)
    Design(nn2, cc2)
end

design_simon = simon(4, 18, 10, 33)
@test abs(expected_sample_size(design_simon, p0) - 22.3) < .1 # precision given in Simons table

obj  = MiniMaxSampleSize(.33, PointMass(p0))
toer = maximal_type_one_error_rate(p0, α)
pwr  = minimal_expected_power(
    PointMass(p1), p0, 1 - β
)

design = Problem(obj, toer, pwr, n1values = collect(10:35), nmax = 35) |>
    pr -> optimise(pr, mle_incompatible_is_error = false)

@test maximum(n(design)) < maximum(n(design_simon))
