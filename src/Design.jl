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
