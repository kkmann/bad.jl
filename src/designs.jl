abstract type AbstractDesign end

# make designs iterable
Base.iterate(design::AbstractDesign, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::AbstractDesign) = 1

Base.show(io::IO, design::AbstractDesign) = print(io, string(design))
Base.show(io::IO, ::MIME"application/prs.juno.inline", design::AbstractDesign) = print(io, string(design))

function string(design::AbstractDesign)
    nn1          = n1(design)
    cont_region = (early_futility(design) + 1):(early_efficacy(design) - 1)
    if length(cont_region) == 0 # one stage design
        return @sprintf "%s<n=%i;c=%i>" string(typeof(design)) nn1 early_futility(design)
    end
    n2_cont     = n2.(design, cont_region)
    if minimum(n2_cont) == maximum(n2_cont) # group sequential design
        return @sprintf "%s<n1=%i;n2(%i...%i)=%i>" string(typeof(design)) nn1 cont_region[1] cont_region[end] n2_cont[1]
    end
    # generic two stage design
    return @sprintf "%s<n1=%i;n2:[%i,%i]->[%i,%i]>" string(typeof(design)) nn1 cont_region[1] cont_region[end] minimum(n2_cont) maximum(n2_cont)
end

n1(design::AbstractDesign)          = length(design.n2) - 1
n2(design::AbstractDesign, x1::Int) = valid(design, x1) ? design.n2[x1 + 1] : error("0 <= x1 <= n1 violated")
n(design::AbstractDesign, x1::Int)  = n1(design) + n2(design, x1)

early_futility(design::AbstractDesign) = any(design.c2 .== EarlyFutility) ? findlast(design.c2 .== EarlyFutility) - 1 : -1
early_efficacy(design::AbstractDesign) = any(design.c2 .== EarlyEfficacy) ? findfirst(design.c2 .== EarlyEfficacy) - 1 : n1(design) + 1
c2(design::AbstractDesign, x1::Int) = valid(design, x1) ? design.c2[x1 + 1] : error("0 <= x1 <= n1 violated")

function probability(x1::Int, x2::Int, design::AbstractDesign, p::T) where {T<:Real}
    !valid(design, x1, x2) ? (return 0.0) : nothing
    return dbinom(x1, n1(design), p)*dbinom(x2, n2(design, x1), p)
end
probability(x1::Int, x2::Int, design::AbstractDesign, Prior::Prior) = expected_value(p -> probability(x1, x2, design, p), prior)

function probability(x1::Int, design::AbstractDesign, p::T) where {T<:Real}
    !valid(design, x1) ? (return 0.0) : nothing
    return dbinom(x1, n1(design), p)
end
probability(x1::Int, design::AbstractDesign, Prior::Prior) = integrate(Prior, probability.(x1, design, Prior.pivots))

function probability(event::Symbol, design::AbstractDesign, p)
    event == :efficacy ? (return power(design, p)) : nothing
    event == :futility ? (return 1 - power(design, p)) : nothing
    x1 = collect(0:n1(design))
    event == :early_efficacy ? (return sum(dbinom(x1 .== Efficacy(), n1(design), p)) ) : nothing
    event == :early_futility ? (return sum(dbinom(x1 .== Futility(), n1(design), p)) ) : nothing
end
probability(event::Symbol, design::AbstractDesign, Prior::Prior) = expected_value(p -> probability(event, design, p), prior)

function power(x1::Int, design::AbstractDesign, p::T) where {T<:Real}
    !valid(design, x1) ? error("invalid x1") : nothing
    power(n2(design, x1), c2(design, x1), p)
end
power(x1::Int, design::AbstractDesign, cprior::Prior) = expected_value(p -> power(x1, design, p), update(prior, x1, n1(design)))

power(design::AbstractDesign, p::T) where {T<:Real} = sum(probability.(0:n1(design), design, p) .* power.(0:n1(design), design, p))
power(design::AbstractDesign, cprior::Prior) = expected_value(p -> power.(design, p), cprior)

sample_space(design::AbstractDesign) =
    [ [x1, x2] for x1 in 0:n1(design), x2 in 0:maximum(design.n2) if valid(design, x1, x2) ] |>
    x -> hcat(x...)' |>
    x -> convert(Array{Int,2}, x)

reject_null(x2::Int, c2::CriticalValue) = x2 > c2

function reject_null(x1::Int, x2::Int, design::AbstractDesign)
    !valid(design, x1, x2) ? error("invalid x1 / x2 for given design") : nothing
    return reject_null(x2, c2(design, x1))
end

valid(design::AbstractDesign, x1::Int) = 0 <= x1 <= n1(design)
valid(design::AbstractDesign, x1::Int, x2::Int) = valid(design, x1) & (0 <= x2 <= n2(design, x1))

function power(p::T, x::Vector{Bool}, design::AbstractDesign) where {T<:Union{Real, Prior}}
    power(p, x, n1(design), design.n2, design.c2)
end

function power(p::T, xx1::Int, nn1::Int, xx2::Int, nn2::Int, design::AbstractDesign) where {T<:Union{Real, Prior}}
    power(p, xx1, nn1, xx2, nn2, n1(design), design.n2, design.c2)
end



mutable struct Design <: AbstractDesign
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



mutable struct OptimalDesign <: AbstractDesign
    n2::Vector{Int}
    c2::Vector{CriticalValue}
    model::DesignIPModel
    score::Real
    # check n2/c2
end
