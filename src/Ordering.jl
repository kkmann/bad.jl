abstract type Ordering end

Base.iterate(ordering::Ordering, state = 0) = state > 0 ? nothing : (ordering, state + 1)
Base.length(ordering::Ordering) = 1

Base.show(io::IO, ordering::Ordering) = print(io, string(ordering))
Base.show(io::IO, ::MIME"application/prs.juno.inline", ordering::Ordering) = print(io, string(ordering))

function p_value(x1, x2, p0, ordering, design)
    space = sample_space(design)
    return sum(probability.(space[:,1], space[:,2], design, p0)[.!ordering.(space[:,1], space[:,2], x1, x2, design)])
end
