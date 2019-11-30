struct PointMass{T<:Real} <: Prior
    atom::T
    PointMass{T}(atom::T) where{T<:Real} = (0 <= atom <= 1) ? new(atom) : error("atom must be in [0, 1]")
end
PointMass(atom::T) where{T<:Real} = PointMass{T}(atom)

function condition(prior::PointMass{T}; low::T = prior.atom, high::T = prior.atom) where {T<:Real}
    (low <= prior.atom <= high) ? (return prior) : error("conditioning only well-defined when probability atom is contained in intervsl")
end

update(prior::PointMass{T}, x::Int, n::Int) where {T<:Real} = !(0 <= x <= n) ? error("invalid x / n") : (return prior)

expected_value(f::Function, prior::PointMass) = f(prior.atom)

mean(prior::PointMass) = prior.atom

predictive_pmf(x, n, prior::PointMass{T}) where {T<:Real} = dbinom.(x, n, prior.atom)

string(prior::PointMass) = @sprintf "PointMass(atom=%.2f)" prior.atom

cdf(prior::PointMass, p::Real) = (p >= prior.atom) ? 1 : 0
