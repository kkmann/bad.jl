"""
    Design{T<:Real}

Binary two-stage design object defining a discrete sampe size function
``n(x_1)`` and ``c(x_1)``, where ``x_1`` is the observed number of stage-one events
and the test rejects the superiority null hypothesis whenever ``x_1 + x_2 > c(x_1)``.
"""
mutable struct Design
    n::Vector{Int}
    c::Vector{CriticalValue}
    function Design(n::Vector{N}, c::Vector{C}) where {N<:Integer, C<:CriticalValue}
        any(n .< (length(n) - 1)) ? error("n must always be larger than n1") : nothing
        length(n) != length(c) ? error("n and c must be of equal length") : nothing
        new(n, c)
    end
end

"""
    Design(n, c)

Construct a [`Design`](@ref) object from `n` and `c` directly.
Here, `n` must be convertable to an integer vector and `c` to a ???
(both of same length ``n_1``).
Then, `n[x1 + 1]` is the final sample size after observing `x1` responses in
stage one and `c[x1 + 1]` the corresponding critical value.

# Parameters

| Parameter  | Description |
| ---------: | :---------- |
| n          | sample size function (must be convertable to integer vector) |
| c          | critical value function (must be convertable to real vector) |
"""
Design(n, c) = Design(convert(Vector{Integer}, n), convert(Vector{CriticalValue}, c))

Base.show(io::IO, design::Design) = print("Design")

# make designs iterable
Base.iterate(design::Design, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Design) = 1

n(design::Design, x1::Int) = valid(design, x1) ? design.n[x1 + 1] : error("0 <= x1 <= n1 violated")
n1(design::Design) = length(design.n) - 1
n2(design::Design, x1::Int) = n(design, x1) - n1(design)

c(design::Design, x1::Int) = valid(design, x1) ? design.c[x1 + 1] : error("0 <= x1 <= n1 violated")



function probability(x1::Int, x2::Int, design::Design, p::T) where {T<:Real}
    if !valid(design, x1)
        return 0.
    else
        nn  = n(design, x1)
        nn2 = n2(design, x1)
        return !(0 <= x2 <= nn2) ? 0 : p^(x1 + x2)*(1 - p)^(nn - x1 - x2)*dbinom(x1, n1(design), p)*dbinom(x2, nn2, p)
    end
end
probability(x1::Int, x2::Int, design::Design, prior::Prior) = integrate(prior, probability.(x1, x2, design, prior.pivots))

function probability(x1::Int, design::Design, p::T) where {T<:Real}
    if !valid(design, x1)
        return 0.
    else
        return dbinom(x1, n1(design), p)
    end
end
probability(x1::Int, design::Design, prior::Prior) = integrate(prior, probability.(x1, design, prior.pivots))

function probability(event::Symbol, design::Design, p)
    if event == :efficacy
        return power(design, p)
    end
    if event == :futility
        return 1 - power(design, p)
    end
    x1 = collect(0:n1(design))
    if event == :early_efficacy
        return sum(dbinom(x1 .== Efficacy(), n1(design), p))
    end
    if event == :early_futility
        return sum(dbinom(x1 .== Futility(), n1(design), p))
    end
end
probability(event::Symbol, design::Design, prior::Prior) = integrate(prior, probability.(event, design, prior.pivots))


function power(x1::Int, design::Design, p::T) where {T<:Real}
    !valid(design, x1) ? error("invalid x1") : nothing
    cc, nn2 = c(design, x1), n2(design, x1)
    if isa(cc, Efficacy)
        return 1.0
    end
    if isa(cc, Futility)
        return 0.0
    end
    return 1 - pbinom(cc - x1, nn2, p)
end
power(design::Design, p::T) where {T<:Real} = sum(probability.(collect(0:n1(design)), design, p) .* power.(collect(0:n1(design)), design, p))

function power(design::Design, prior::Prior; mrv::T = 0.) where {T<:Real}
    cprior = condition(prior, low = mrv)
    return integrate(cprior, power.(design, cprior.pivots) )
end

function power(x1::Int, design::Design, prior::Prior; mrv::T = 0.) where {T<:Real}
    !valid(design, x1) ? error("invalid x1") : nothing
    cposterior = condition(posterior(prior, x1, n1(design)), low = mrv)
    return integrate(cposterior, power.(x1, design, cposterior.pivots) )
end



function reject_null(x1::Int, x2::Int, design::Design)
    !valid(design, x1, x2) ? error("invalid x1 / x2for given design") : nothing
    return x1 + x2 > c(design, x1)
end


# """
#     simulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}
# Simulate `nsim` trials of design `design` with true response probability `p`.
# # Parameters
# | Parameter    | Description |
# | -----------: | :---------- |
# | design       | a [`Design`](@ref) |
# | p            | true response probability |
# | nsim         | number of simulated trials |
# # Return Value
# A `DataFrame` with columns `x1` (stage one responses) `n` (overall sample size)
# `c` (critical value), `x2` (stage-two responses), and `rejectedH0` (test decision).
# """
# function simulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}
#
#   x2    = SharedArray{Int}(nsim)
#   n     = SharedArray{Int}(nsim)
#   c     = SharedArray{Float64}(nsim)
#   rej   = SharedArray{Bool}(nsim)
#   n1    = interimsamplesize(design)
#   rv_x1 = Distributions.Binomial(n1, p)
#   x1    = SharedArray{Int}(nsim)
#   @sync @parallel for i in 1:nsim
#       x1[i]  = rand(rv_x1)
#       n2     = samplesize(design, x1[i]) - n1
#       n[i]   = n2 + n1
#       c[i]   = criticalvalue(design, x1[i])
#       x2[i]  = rand(Distributions.Binomial(n2, p))
#       rej[i] = test(design, x1[i], x2[i])
#   end
#   return DataFrames.DataFrame(
#       x1 = convert(Vector{Int}, x1),
#       n  = convert(Vector{Int}, n),
#       c  = convert(Vector{Float64}, c),
#       x2 = convert(Vector{Int}, x2),
#       rejectedH0 = convert(Vector{Bool}, rej)
#   )
#
# end

# """
#     writecsv(filename::String, design::Design; label::String = "")
# Save design as .csv file.
# # Parameters
# | Parameter    | Description |
# | -----------: | :---------- |
# | filename     | filename for sving design |
# | design       | a [`Design`](@ref) |
# | label        | string name of design |
# """
# function writecsv(filename::String, design::Design; label::String = "")
#
#   df = DataFrames.DataFrame(design)
#   label != "" ? df[:design] = label : nothing
#   CSV.write(filename, df)
#   return df
#
# end


function get_x1_x2_grid(design::Design)
    nmax = maximum(design.n)
    return [[x1, x2] for x1 in 0:n1(design), x2 in 0:(nmax - n1(design)) if
        (x2 <= n(design, x1) - n1(design)) & ((n(design, x1) > n1(design)) | (x2 == 0))
    ] |> x -> hcat(x...)'
end

valid(design::Design, x1::Int) = 0 <= x1 <= n1(design)
valid(design::Design, x1::Int, x2::Int) = valid(design, x1) & (0 <= x2 <= n(design, x1) - n1(design))
