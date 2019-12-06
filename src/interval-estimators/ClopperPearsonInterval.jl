struct ClopperPearsonInterval{TO<:Ordering,TD<:AbstractDesign,TI<:Integer,TR<:Real} <: IntervalEstimator
    ordering::TO
    design::TD
    XX::Array{TI,2}
    bounds::Array{TR,2}
    alpha::TR
end

string(ci::ClopperPearsonInterval) = @sprintf "ClopperPearsonInterval<%s>" string(ci.ordering)

function ClopperPearsonInterval(ordering::Ordering, design::AbstractDesign, α::Real; ϵ = 1e-6)

    XX = sample_space(design)
    nn = size(XX, 1)

    pval_sup(p0, x1, x2) = p_value(x1, x2, p0, ordering, design; orientation = :superiority)
    pval_inf(p0, x1, x2) = p_value(x1, x2, p0, ordering, design; orientation = :inferiority)

    bounds = zeros(Float64, nn, 2)
    for i in 1:nn
        x1, x2 = XX[i,:]
        phat   = (x1 + x2) / n(design, x1)
        l = try
                Roots.find_zero(p -> pval_sup(p, x1, x2) - (α + ϵ), phat)
            catch e
                isa(e, Roots.ConvergenceFailed) ? 0.0 : e
        end
        u = try
                Roots.find_zero(p -> pval_inf(p, x1, x2) - (α + ϵ), phat)
            catch e
                isa(e, Roots.ConvergenceFailed) ? 1.0 : e
        end
        u <= l ? error("confidence interval collapsed") : nothing
        bounds[i,:] .= [l, u]
    end
    return ClopperPearsonInterval{typeof(ordering),typeof(design),eltype(XX),Float64}(ordering, design, XX, bounds, α)
end
