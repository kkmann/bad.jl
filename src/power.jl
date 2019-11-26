prune(p) = min(1, max(0, p))

power(n::Int, c::Real, p::T) where {T<:Real} =
    1 - pbinom(c, n, p)
power(n::Int, c::Real, prior::Prior) =
    expected_value(p -> power(n, c, p), prior) |> prune
power(x1::Int, n1::Int, n2::Int, c2::Real, prior::Prior) =
    expected_value(p -> power(n2, c2, p), update(prior, x1, n1)) |> prune


function power(x1::Int, design::AbstractDesign, p::T) where {T<:Real}
    !valid(design, x1) ? error("invalid x1") : nothing
    power(n2(design, x1), c2(design, x1), p) |> prune
end
power(x1::Int, design::AbstractDesign, cprior::Prior) =
    expected_value(p -> power(x1, design, p), update(cprior, x1, n1(design))) |> prune

power(design::AbstractDesign, p::T) where {T<:Real} =
    sum(dbinom.(0:n1(design), n1(design), p) .* power.(0:n1(design), design, p)) |> prune
power(design::AbstractDesign, cprior::Prior) =
    expected_value(p -> power(design, p), cprior) |> prune



function power(p::Real, xx1::Int, nn1::Int, n1::Int, n2::Vector{Int}, c2::Vector{T}) where {T<:Real}
    any(length.((n2, c2)) .!= n1 + 1) ? error("n2/c2 must be of length n1 + 1") : nothing
    nn1 > n1 ? error("") : nothing
    x1_delta = 0:(n1 - nn1)
    x1       = xx1 .+ x1_delta
    return sum( dbinom.(x1_delta, n1 - nn1, p) .* power.(n2[x1 .+ 1], c2[x1 .+ 1], p) ) |> prune
end
function power(p::Prior, xx1::Int, nn1::Int, n1::Int, n2::Vector{Int}, c2::Vector{T}) where {T<:Real}
    any(length.((n2, c2)) .!= n1 + 1) ? error("n2/c2 must be of length n1 + 1") : nothing
    nn1 > n1 ? error("") : nothing
    p = update(p, xx1, nn1)
    x1_delta = 0:(n1 - nn1)
    x1       = xx1 .+ x1_delta
    return sum( dbinom.(x1_delta, n1 - nn1, p) .* power.(n2[x1 .+ 1], c2[x1 .+ 1], p) ) |> prune
end
function power(p::T, xx1::Int, nn1::Int, design::AbstractDesign) where {T<:Union{Real,Prior}}
    return power(p, xx1, nn1, n1(design), design.n2, design.c2) |> prune
end
