struct MarginalFeasibleSpace{TI<:Integer,TR<:Real}
    n1values::Vector{TI}
    nmax::TI
    n2mincontinueabs::TI
    n2mincontinuereln1::TR
    n2maxcontinuereln1::TR
    type::Symbol
    x1min::TI
    pnull::TR
    α::TR
    palt::TR
    β::TR
end
function MarginalFeasibleSpace(
        pnull, palt;
        maxmultipleonestage = 2.,
        α = .025, β = .1,
        nmax = Int(ceil(maxmultipleonestage*one_stage_sample_size(pnull, α, palt, β))),
        n1minfctr = .25,
        n1min = Int(ceil(n1minfctr*one_stage_sample_size(pnull, α, palt, β))),
        n1maxfctr = .66,
        n1max = Int(ceil(n1maxfctr*nmax)),
        n1values = collect(nmin:n1max),
        n2mincontinueabs = 5,
        n2mincontinuereln1 = 1.1,
        n2maxcontinuereln1 = 4.,
        type = :TwoStage,
        x1min = 0
    )
    return MarginalFeasibleSpace{typeof(nmax),typeof(pnull)}(
        n1values, nmax,
        n2mincontinueabs, n2mincontinuereln1, n2maxcontinuereln1,
        type,
        x1min,
        pnull,
        α,
        palt,
        β
    )
end

n1(mfs::MarginalFeasibleSpace) = mfs.n1values
