abstract type Prior end

# make priors iterable
Base.iterate(design::Prior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Prior) = 1

Base.show(io::IO, prior::Prior) = print(io, string(prior))
Base.show(io::IO, ::MIME"application/prs.juno.inline", prior::Prior) = print(io, string(prior))
