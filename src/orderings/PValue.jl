struct PValue{TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign}
    ordering::TO
    design::TD
    XX::Array{TI,2}
    inds::Array{Bool,2}
    p0::TR
end

Base.iterate(pval::PValue, state = 0) = state > 0 ? nothing : (pval, state + 1)
Base.length(pval::PValue) = 1

Base.show(io::IO, pval::PValue) = print(io, string(pval))
Base.show(io::IO, ::MIME"application/prs.juno.inline", pval::PValue) = print(io, string(pval))
string(pval::PValue) = "p-value"



function PValue(ordering::Ordering, design::AbstractDesign, p0::Real)

    XX   = sample_space(design)
    nn   = size(XX, 1)
    inds = zeros(nn, nn)
    for i = 1:nn
        x1, x2    = XX[i,1], XX[i,2]
        inds[i,:] = more_extreme.(XX[:,1], XX[:,2], x1, x2, ordering, design)
    end
    return PValue{eltype(XX),typeof(p0),typeof(ordering),typeof(design)}(ordering, design, XX, inds, p0)
end

function PValue(estimator::Estimator, design::AbstractDesign, p0::Real; orientation::Symbol = :superiority)
    PValue(EstimatorOrdering(estimator), design, p0)
end



function (pval::PValue{TI,TR,TO,TD})(x1::TI, x2::TI) where {TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign}

    design = pval.design
    p0     = pval.p0
    XX     = pval.XX
    ind    = mapslices(all, XX .== [x1 x2], dims = 2) |> vec |> findall
    @assert length(ind) == 1
    ind    = ind[1]
    inds   = pval.inds[ind,:]
    return min(1, max(0, sum(
        pmf.(XX[inds,2], n2.(design, XX[inds,1]), p0) .*
        pmf.(XX[inds,1], n1(design), p0)
    ) ) )
end

evaluate(pval::PValue{TI,TR,TO,TD}, x1::TI, x2::TI) where {TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign} = pval(x1, x2)
