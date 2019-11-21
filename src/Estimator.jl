abstract type Estimator end

Base.iterate(design::Estimator, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Estimator) = 1

Base.show(io::IO, design::Estimator) = print(io, string(design))
Base.show(io::IO, ::MIME"application/prs.juno.inline", design::Estimator) = print(io, string(design))
