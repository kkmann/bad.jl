abstract type AbstractDesign end

# make designs iterable
Base.iterate(design::AbstractDesign, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::AbstractDesign) = 1

Base.show(io::IO, ::MIME"application/prs.juno.inline", design::AbstractDesign) = print(io, string(design))

function string(design::AbstractDesign)
    nn1          = n1(design)
    cont_region = (early_futility(design) - 1):(early_efficacy(design) + 1)
    if length(cont_region) == 0 # one stage design
        return @sprintf "%s<n=%i;c=%i>" string(typeof(design)) nn1 early_futility(design)
    end
    n2_cont     = n2.(design, cont_region)
    if minimum(n2_cont) == maximum(n2_cont) # group sequential design
        return @sprintf "%s<n1=%i;n2:[%i,%i]=%i>" string(typeof(design)) nn1 cont_region[1] cont_region[end] n2_cont[1]
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
    !valid(design, x1) ? (return 0.0) : nothing
    nn, nn2 = n(design, x1), n2(design, x1)
    return !(0 <= x2 <= nn2) ? 0 : p^(x1 + x2)*(1 - p)^(nn - x1 - x2) * dbinom(x1, n1(design), p) * dbinom(x2, nn2, p)
end
probability(x1::Int, x2::Int, design::AbstractDesign, Prior::Prior) = integrate(Prior, probability.(x1, x2, design, Prior.pivots))

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
probability(event::Symbol, design::AbstractDesign, Prior::Prior) = integrate(Prior, probability.(event, design, Prior.pivots))

function power(x1::Int, n2::Int, c2::CriticalValue, p::T) where {T<:Real} # does not depend on x1
    isa(c2, Efficacy) ? (return 1.0) : nothing
    isa(c2, Futility) ? (return 0.0) : nothing
    return 1 - pbinom(c2, n2, p)
end
power(x1::Int, n1::Int, n2::Int, c2::CriticalValue, cprior::Prior) = integrate(update(cprior, x1, n1), power.(x1, n2, c2, cprior.pivots))

function power(x1::Int, design::AbstractDesign, p::T) where {T<:Real}
    !valid(design, x1) ? error("invalid x1") : nothing
    power(x1, n2(design, x1), c2(design, x1), p)
end
power(x1::Int, design::AbstractDesign, cprior::Prior) = integrate(update(cprior, x1, n1(design)), power.(x1, design, cprior.pivots) )

power(design::AbstractDesign, p::T) where {T<:Real} = sum(probability.(0:n1(design), design, p) .* power.(0:n1(design), design, p))
power(design::AbstractDesign, cprior::Prior) = integrate(cprior, power.(design, cprior.pivots) )



function reject_null(x1::Int, x2::Int, design::AbstractDesign)
    !valid(design, x1, x2) ? error("invalid x1 / x2 for given design") : nothing
    return x2 > c2(design, x1)
end

valid(design::AbstractDesign, x1::Int) = 0 <= x1 <= n1(design)
valid(design::AbstractDesign, x1::Int, x2::Int) = valid(design, x1) & (0 <= x2 <= n2(design, x1))
