using Test, bad; import Plots

# get curtailed simons design (implied stopping for efficacy!)
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

function optimal_design(pnull, palt)
    essnull = SampleSize(p | pnull)
    problem = Problem(
        minimise(essnull),
        subject_to(TypeOneErrorRate(p | pnull), α, [-.1, 1.1]),
        subject_to(Power(p | palt), β, [-.1, 1.1]),
        nmax                    = convert(Int, ceil(1.5*maximum(n(simonsdesign)))),
        n1values                = convert(Vector{Int}, floor(.5*n1(simonsdesign)):ceil(1.5*n1(simonsdesign))),
        # relax all other heuristics
        n2maxcontinuereln1      = 10.,
        n2mincontinuereln1      = 0.0,
        n2mincontinueabs        = 5,
        curtail_stage_one_fct   = 0.0,
        type                    = :GroupSequential
    )
    design = optimise(problem; verbosity = 0)
    return design, problem, size(problem)
end


α, β  = 0.05, .2
p     = Beta(1, 1) # uniform prior, only used to make conditioning valid for all p ∈ [0, 1]

# group sequential designs are notoriously hard to fit due to the global
# constraint, only consider the two cornering cases (low/high p0)
# and one in the middle, values taken from:
#
# > Simon, R. (1989). Optimal two-stage designs for phase II clinical trials.
# > Controlled clinical trials, 10(1), 1-10.

pnull, palt  = .05, .25; simonsdesign = simon(0, 9, 2, 17)
design, problem, nvar = optimal_design(pnull, palt)
@test all( problem.toer.score.((design, simonsdesign)) .<= α )
@test all( problem.power.score.((design, simonsdesign)) .>= 1 - β )
@test problem.objective(design) <= problem.objective(simonsdesign)


pnull, palt  = .1, .3; simonsdesign = simon(1, 10, 5, 29)
design, problem, nvar = optimal_design(pnull, palt)
@test all( problem.toer.score.((design, simonsdesign)) .<= α )
@test all( problem.power.score.((design, simonsdesign)) .>= 1 - β )
@test problem.objective(design) <= problem.objective(simonsdesign)


pnull, palt  = .7, .9; simonsdesign = simon(4, 6, 22, 27)
design, problem, nvar = optimal_design(pnull, palt)
@test all( problem.toer.score.((design, simonsdesign)) .<= α )
@test all( problem.power.score.((design, simonsdesign)) .>= 1 - β )
@test problem.objective(design) <= problem.objective(simonsdesign)
