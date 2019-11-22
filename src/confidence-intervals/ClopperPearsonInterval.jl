struct ClopperPearsonInterval{T<:Real} <: ConfidenceInterval
    p::Vector{T}
    p_values_superiority::Matrix{T}
    p_values_inferiority::Matrix{T}
    space::Matrix{Int}
    ordering::Ordering
    design::AbstractDesign
end

string(ci::ClopperPearsonInterval) = @sprintf "ClopperPearsonInterval<%s>" string(ci.ordering)

function ClopperPearsonInterval(ordering::Ordering, design::AbstractDesign; resolution = .001)
    p                    = collect(0:.001:1)
    space                = sample_space(design)
    p_values_superiority = zeros(size(space, 1), length(p))
    p_values_inferiority = zeros(size(space, 1), length(p))
    for i in 1:length(p)
        p_values_superiority[:, i] = p_value.(space[:, 1], space[:, 2], p[i], ordering, design; orientation = :superiority)
        p_values_inferiority[:, i] = p_value.(space[:, 1], space[:, 2], p[i], ordering, design; orientation = :inferiority)
    end
    ClopperPearsonInterval(p, p_values_superiority, p_values_inferiority, space, ordering, design)
end

function (ci::ClopperPearsonInterval{T})(x1::Int, x2::Int; confidence::T = .9) where{T<:Real}
    threshold = (1 - confidence)/2
    # find index in sample space
    inds = [ [x1, x2] == ci.space[i,:] for i in 1:size(ci.space, 1) ]
    sum(inds) == 0 ? error("x1/x2 not found in sample space, valid? ") : nothing
    sum(inds)  > 1 ? error("x1/x2 not unique in sample space, check sample space? ") : nothing
    ind = findfirst(inds)
    hi = maximum(ci.p[ci.p_values_inferiority[ind, :] .> threshold])
    lo = minimum(ci.p[ci.p_values_superiority[ind, :] .> threshold])
    lo > hi ? error("something went wrong") : nothing
    return [lo, hi]
end
