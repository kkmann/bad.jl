struct PValue{TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign}
    ordering::TO
    design::TD
    ordered_sample_space::Array{TI,2}
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
    rank = zeros(nn)
    for i = 1:nn
        x1, x2  = XX[i,1], XX[i,2]
        rank[i] = sum( .!more_extreme.(XX[:,1], XX[:,2], x1, x2, ordering, design) )
    end
    XX = XX[sortperm(rank), :]
    return PValue{eltype(XX),typeof(p0),typeof(ordering),typeof(design)}(ordering, design, XX, p0)
end

function PValue(estimator::Estimator, design::AbstractDesign, p0::Real; orientation::Symbol = :superiority)
    PValue(EstimatorOrdering(estimator), design, p0)
end



function (pval::PValue{TI,TR,TO,TD})(x1::TI, x2::TI) where {TI<:Integer,TR<:Real,TO<:Ordering,TD<:AbstractDesign}
    for i in 1:size(pval.ordered_sample_space, 1)
        if [x1, x2] == pval.ordered_sample_space[i,:]
            return min(1, max(0, sum(
                pdf.(
                    pval.ordered_sample_space[i:end, 1],
                    pval.ordered_sample_space[i:end,2],
                    pval.design,
                    pval.p0
                )
            ) ) )
        end
    end
    error("(x1,x2) not found in sample space, valid observation?")
end
