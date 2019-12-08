struct MarginalFeasibleSpace{TI<:Integer,TR<:Real}
    n1values::Vector{TI}
    nmax::TI
    n2mincontinueabs::TI
    n2mincontinuereln1::TR
    n2maxcontinuereln1::TR
    type::Symbol
    x1min::TI
end
function MarginalFeasibleSpace(
        n1values, nmax;
        n2mincontinueabs = 5, n2mincontinuereln1 = 1.1, n2maxcontinuereln1 = 4.,
        type = :TwoStage, x1min = 0
    )
    return MarginalFeasibleSpace{typeof(nmax),typeof(n2mincontinuereln1)}(
        n1values, nmax;
        n2mincontinueabs   = n2mincontinueabs,
        n2mincontinuereln1 = n2mincontinuereln1,
        n2maxcontinuereln1 = n2maxcontinuereln1,
        type               = type,
        x1min              = x1min
    )
)
