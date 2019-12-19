using Test, bad; import Plots


α, β        = 0.05, .2
pnull, palt = .2, .4
n_required  = 35 # for corresponding single stage test
p           = Beta(1, 1) # uniform prior, only used to make conditioning valid for all p ∈ [0, 1]

# check under expected sample sise
problem = Problem(
    minimise(
        SampleSize(p | palt)
    ),
    Power(p | pnull) <= α,
    Power(p | palt)  >= 1 - β,
    # make search space large enough to incorporate single-stage solution
    nmax                    = n_required + 10,
    n1max                   = n_required + 1,
    type                    = :OneStage
)
design = optimise(problem; verbosity = 0)
@test problem.toer.score(design) .<= α
@test problem.power.score(design) .>= 1 - β
@test all( n2.(design, 0:n1(design)) .== 0 )
@test n1(design) == n_required


# check under minimax sample size
problem = Problem(
    MiniMaxSampleSize(.33, p),
    Power(p | pnull) <= α,
    Power(p | palt)  >= 1 - β,
    # make search space large enough to incorporate single-stage solution
    nmax                    = n_required + 10,
    n1max                   = n_required + 1,
    type                    = :OneStage
)
design = optimise(problem; verbosity = 0)
@test problem.toer.score(design) .<= α
@test problem.power.score(design) .>= 1 - β
@test all( n2.(design, 0:n1(design)) .== 0 )
@test n1(design) == n_required
