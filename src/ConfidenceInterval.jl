abstract type ConfidenceInterval end

Base.iterate(ci::ConfidenceInterval, state = 0) = state > 0 ? nothing : (ci, state + 1)
Base.length(ci::ConfidenceInterval) = 1

Base.show(io::IO, ci::ConfidenceInterval) = print(io, string(ci))
Base.show(io::IO, ::MIME"application/prs.juno.inline", ci::ConfidenceInterval) = print(io, string(ci))
