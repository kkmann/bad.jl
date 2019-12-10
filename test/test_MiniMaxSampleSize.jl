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
    # cc2[cc2 .< 0]     .= -Inf
    # cc2[cc2 .>= nn2]  .= Inf
    # nn2[cc2 .== -Inf] .= 0
    Design(nn2, cc2)
end

design_simon = simon(4, 18, 10, 33)
ess = SampleSize(PointMass(p0))
@test abs(ess(design_simon) - 22.3) < .1 # precision given in Simons table



minimax = MiniMaxSampleSize(.2, PointMass(p0))
toer    = TypeOneErrorRate(PointMass(p0), p0)
pow     = Power(PointMass(p1), p0 + .1)


problem = Problem(
        minimax,
        subject_to(toer, α, [-.1, .5]),
        subject_to(pow, β, [.5, 1.1]);
        n1values = collect(10:35),
        nmax = 35,
        n2mincontinueabs = 5,
        n2mincontinuereln1 = .9,
        n2maxcontinuereln1 = 5.,
        curtail_stage_one_fct = .01
    )
size(problem)
design = optimise(problem, timelimit = 300)

@test maximum(n(design)) < maximum(n(design_simon))

plot(design)
