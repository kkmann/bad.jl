abstract type Ordering end

Base.iterate(ordering::Ordering, state = 0) = state > 0 ? nothing : (ordering, state + 1)
Base.length(ordering::Ordering) = 1

Base.show(io::IO, ordering::Ordering) = print(io, string(ordering))
Base.show(io::IO, ::MIME"application/prs.juno.inline", ordering::Ordering) = print(io, string(ordering))
