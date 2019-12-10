abstract type Prior end

# make priors iterable
Base.iterate(design::Prior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Prior) = 1

Base.show(io::IO, prior::Prior) = print(io, string(prior))
Base.show(io::IO, ::MIME"application/prs.juno.inline", prior::Prior) = print(io, string(prior))

bounds(prior::Prior) = [prior.low, prior.high]

pdf(p::Real, prior::Prior) = error("not implemented")

cdf(p::Real, prior::Prior) = error("not implemented")

expectation(f::Function, prior::Prior) = error("not implemented")

mean(prior::Prior) = error("not implemented")

condition(prior::Prior; low::Real = prior.low, high::Real = prior.high) = error("not implemented")
function |(prior::Prior, interval::Vector{T}) where {T<:Real}

    @assert (length(interval) == 2) @sprintf "must condition on interval (vector of length 2, is %i)" length(interval)
    @assert (0 <= interval[1] < interval[2] <= 1)  @sprintf "invalid interval specification must be 0<=a<b<=1, is [%.2f, %.2f]" interval[1] interval[2]
    condition(prior; low = interval[1], high = interval[2])
end
condition(prior::Prior, point::Real) = (pdf(point, prior) > 0) ? PointMass(point) : error("conditioning on point with pdf = 0 is not well-defined")
|(prior::Prior, point::Real) = condition(prior, point)

<=(low::T, prior::Prior) where {T<:Real} = condition(prior; low = low)
<=(prior::Prior, high::T) where {T<:Real} = condition(prior; high = high)
>=(prior::Prior, low::T) where {T<:Real} = condition(prior; low = low)
>=(high::T, prior::Prior) where {T<:Real} = condition(prior; high = high)

update(prior::Prior, x::Integer, n::Integer) = error("not implemented")

string(prior::Prior) = error("not implemented")
