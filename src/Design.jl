mutable struct Design
    n2::Vector{Int}
    c2::Vector{CriticalValue}
    function Design(n2::Vector{N2}, c2::Vector{C2}) where {N2<:Integer, C2<:CriticalValue}
        !all(0 .<= n2) ? error("n must be positive") : nothing
        # ToDo: check for ineffective stopping (n2 > 0 and early stopping)
        length(n2) != length(c2) ? error("n and c must be of equal length") : nothing
        new(n2, c2)
    end
end
Design(n2, c2) = Design(convert(Vector{Integer}, n2), convert(Vector{CriticalValue}, c2))

Base.show(io::IO, design::Design) = print("Design")

# make designs iterable
Base.iterate(design::Design, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Design) = 1

n1(design::Design)          = length(design.n2) - 1
n2(design::Design, x1::Int) = valid(design, x1) ? design.n2[x1 + 1] : error("0 <= x1 <= n1 violated")
n(design::Design, x1::Int)  = n1(design) + n2(design, x1)

c2(design::Design, x1::Int) = valid(design, x1) ? design.c2[x1 + 1] : error("0 <= x1 <= n1 violated")

function probability(x1::Int, x2::Int, design::Design, p::T) where {T<:Real}
    !valid(design, x1) ? (return 0.0) : nothing
    nn, nn2 = n(design, x1), n2(design, x1)
    return !(0 <= x2 <= nn2) ? 0 : p^(x1 + x2)*(1 - p)^(nn - x1 - x2) * dbinom(x1, n1(design), p) * dbinom(x2, nn2, p)
end
probability(x1::Int, x2::Int, design::Design, Prior::Prior) = integrate(Prior, probability.(x1, x2, design, Prior.pivots))

function probability(x1::Int, design::Design, p::T) where {T<:Real}
    !valid(design, x1) ? (return 0.0) : nothing
    return dbinom(x1, n1(design), p)
end
probability(x1::Int, design::Design, Prior::Prior) = integrate(Prior, probability.(x1, design, Prior.pivots))

function probability(event::Symbol, design::Design, p)
    event == :efficacy ? (return power(design, p)) : nothing
    event == :futility ? (return 1 - power(design, p)) : nothing
    x1 = collect(0:n1(design))
    event == :early_efficacy ? (return sum(dbinom(x1 .== Efficacy(), n1(design), p)) ) : nothing
    event == :early_futility ? (return sum(dbinom(x1 .== Futility(), n1(design), p)) ) : nothing
end
probability(event::Symbol, design::Design, Prior::Prior) = integrate(Prior, probability.(event, design, Prior.pivots))

function power(x1::Int, n2::Int, c2::CriticalValue, p::T) where {T<:Real} # does not depend on x1
    isa(c2, Efficacy) ? (return 1.0) : nothing
    isa(c2, Futility) ? (return 0.0) : nothing
    return 1 - pbinom(c2, n2, p)
end
power(x1::Int, n1::Int, n2::Int, c2::CriticalValue, cprior::Prior) = integrate(update(cprior, x1, n1), power.(x1, n2, c2, cprior.pivots))

function power(x1::Int, design::Design, p::T) where {T<:Real}
    !valid(design, x1) ? error("invalid x1") : nothing
    power(x1, n2(design, x1), c2(design, x1), p)
end
power(x1::Int, design::Design, cprior::Prior) = integrate(update(cprior, x1, n1(design)), power.(x1, design, cprior.pivots) )

power(design::Design, p::T) where {T<:Real} = sum(probability.(0:n1(design), design, p) .* power.(0:n1(design), design, p))
power(design::Design, cprior::Prior) = integrate(cprior, power.(design, cprior.pivots) )



function reject_null(x1::Int, x2::Int, design::Design)
    !valid(design, x1, x2) ? error("invalid x1 / x2for given design") : nothing
    return x2 > c2(design, x1)
end

valid(design::Design, x1::Int) = 0 <= x1 <= n1(design)
valid(design::Design, x1::Int, x2::Int) = valid(design, x1) & (0 <= x2 <= n2(design, x1))
