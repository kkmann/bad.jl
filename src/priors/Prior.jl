abstract type Prior end

# make priors iterable
Base.iterate(design::Prior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Prior) = 1

Base.show(io::IO, prior::Prior) = print(io, string(prior))
Base.show(io::IO, ::MIME"application/prs.juno.inline", prior::Prior) = print(io, string(prior))

pdf(p::Real, prior::Prior) = error("not implemented")

# pdf(x1::Integer, x2::Integer, design::AbstractDesign, prior::Prior) = error("not implemented")

cdf(p::Real, prior::Prior) = error("not implemented")

expectation(f::Function, prior::Prior) = error("not implemented")

mean(prior::Prior) = error("not implemented")

condition(prior::Prior; low::Real = prior.low, high::Real = prior.high) = error("not implemented")

update(prior::Prior, x::Integer, n::Integer) = error("not implemented")

string(prior::Prior) = error("not implemented")
