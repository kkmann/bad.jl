struct PointMass{T<:Real} <: Prior
    atom::T
    PointMass{T}(atom::T) where{T<:Real} = (0 <= atom <= 1) ? new(atom) : error("atom must be in [0, 1]")
end
PointMass(atom::T) where{T<:Real} = PointMass{T}(atom)

pdf(prior::PointMass{T}, p::T) where{T<:Real} = (p == prior.atom) ? Inf : 0

cdf(prior::PointMass{T}, p::T) where{T<:Real} = (p >= prior.atom) ? 1 : 0

function condition(prior::PointMass{T}; low::T = prior.atom, high::T = prior.atom) where {T<:Real}
    (low <= prior.atom <= high) ? (return prior) : error("conditioning only well-defined when probability atom is contained in intervsl")
end

update(prior::PointMass{T}, x::TI, n::TI) where {T<:Real,TI<:Integer} = !(0 <= x <= n) ? error("invalid x / n") : (return prior)

expectation(f::Function, prior::PointMass{T}) where {T<:Real} = f(prior.atom)

mean(prior::PointMass{T}) where {T<:Real} = prior.atom

string(prior::PointMass) = @sprintf "PointMass(%.2f)" prior.atom

quantile(prior::PointMass, prob::Real) = prior.atom
