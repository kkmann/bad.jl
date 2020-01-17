using Test, bad



# Simons minimax design
# (Optimal Two-Stage Designs for Phase II Clinical Trials, Richard Simon)
# r1/n1 = 4/18; r/n = 10/33
function simon(r1::Int, n1::Int, r::Int, n::Int)
    nn2 = zeros(Int, n1 + 1)
    nn2[(r1 + 2):end] .= n - n1
    cc = repeat([Inf], n1 + 1)
    cc[(r1 + 2):end] .= r
    cc2 = cc .- collect(0:n1)
    Design(nn2, cc2)
end



pnull, palt  = .2, .4
α, β         = .05, .2
p            = Beta(1, 1)

design_simon = simon(4, 18, 10, 33)
essnull      = SampleSize(p | pnull)
@test essnull(design_simon) ≈ 22.3 atol = .1 # precision given in Simons table

problem = Problem(
        MiniMaxSampleSize(.2, p | pnull),
        Power(p | pnull) <= α,
        Power(p | palt)  >= 1 - β;
        # again, just removing as many as possible heuristics from the marginal fesible space
        n1values = collect(10:35),
        nmax = 35,
        n2mincontinueabs = 5,
        n2ton1fctrs      = (.9, 5.0),
        curtail_stage_one_buffer = 5
    )
design = optimise(problem; verbosity = 0)

# show that its actually better than the group sequential one ;)
@test maximum(n(design)) < maximum(n(design_simon))
