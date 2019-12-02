abstract type AbstractDesign end

# make designs iterable
Base.iterate(design::AbstractDesign, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::AbstractDesign) = 1

Base.show(io::IO, design::AbstractDesign) = print(io, string(design))
Base.show(io::IO, ::MIME"application/prs.juno.inline", design::AbstractDesign) = print(io, string(design))

function string(design::AbstractDesign)
    nn1         = n1(design)
    cont_region = Int(max(0, early_futility(design) + 1)):Int(min(nn1, early_efficacy(design) - 1))
    if length(cont_region) == 0 # one stage design
        return @sprintf "%s<n=%i;c=%i>" string(typeof(design)) nn1 early_futility(design)
    end
    n2_cont     = n2.(design, cont_region)
    if minimum(n2_cont) == maximum(n2_cont) # group sequential design
        return @sprintf "%s<n1=%i;n2:[%i,%i]->%i>" string(typeof(design)) nn1 cont_region[1] cont_region[end] n2_cont[1]
    end
    # generic two stage design
    return @sprintf "%s<n1=%i;n2:[%i,%i]->[%i,%i]>" string(typeof(design)) nn1 cont_region[1] cont_region[end] minimum(n2_cont) maximum(n2_cont)
end

n1(design::AbstractDesign)             = length(design.n2) - 1
n2(design::AbstractDesign, x1::Int)    = valid(design, x1) ? design.n2[x1 + 1] : error("0 <= x1 <= n1 violated")
n(design::AbstractDesign, x1::Int)     = n1(design) + n2(design, x1)
n(design::AbstractDesign)              = n1(design) .+ design.n2

early_futility(design::AbstractDesign) = any(design.c2 .== Inf) ? findlast(design.c2 .== Inf) - 1 : -Inf
early_efficacy(design::AbstractDesign) = any(design.c2 .== -Inf) ? findfirst(design.c2 .== -Inf) - 1 : Inf
c2(design::AbstractDesign, x1::Int)    = valid(design, x1) ? design.c2[x1 + 1] : error("0 <= x1 <= n1 violated")

expected_sample_size(design::AbstractDesign, p::Real) = sum(dbinom.(0:n1(design), n1(design), p) .* design.n2) + n1(design)

function pdf(x1::TI, x2::TI, design::TD, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign}

    return dbinom(x1, n1(design), p) * dbinom(x2, n2(design,x1), p)
end

function as_table(design::AbstractDesign)
    DataFrames.DataFrame(
        n1 = repeat([n1(design)], n1(design) + 1),
        x1 = 0:n1(design),
      phat = (0:n1(design)) ./ repeat([n1(design)], n1(design) + 1),
        n2 = design.n2,
        n  = n(design),
        c2 = design.c2
    )
end

function plot(design::AbstractDesign)
    tbl = as_table(design)
    Gadfly.plot(
        Gadfly.layer(tbl, x = :phat, y = :n, Gadfly.Geom.hair, Gadfly.Geom.point),
        Gadfly.layer(tbl, x = :phat, y = tbl[!, :c2] + tbl[!, :x1], Gadfly.Geom.hair, Gadfly.Geom.point, Gadfly.Theme(default_color = Gadfly.colorant"orange"))
    )
end



sample_space(design::AbstractDesign) =
    [ [x1, x2] for x1 in 0:n1(design), x2 in 0:maximum(design.n2) if valid(design, x1, x2) ] |>
    x -> hcat(x...)' |>
    x -> convert(Array{Int,2}, x)

reject_null(x2::Int, c2::Real) = x2 > c2

function reject_null(x1::Int, x2::Int, design::AbstractDesign)
    !valid(design, x1, x2) ? error("invalid x1 / x2 for given design") : nothing
    return reject_null(x2, c2(design, x1))
end

valid(design::AbstractDesign, x1::Int) = 0 <= x1 <= n1(design)
valid(design::AbstractDesign, x1::Int, x2::Int) = valid(design, x1) & (0 <= x2 <= n2(design, x1))


function fisher_information(p::T, design::TD)::T where {T<:Real,TD<:AbstractDesign}
    # 0 < p < 1 ? nothing : error("p must be in (0, 1)")
    (p ≈ 0) | (p ≈ 1) ? (return Inf) : nothing
    x1_x2     = sample_space(design)
    log_prob  = log.(dbinom.(x1_x2[:, 1], n1(design), p)) .+ log.(dbinom.(x1_x2[:, 2], n2.(design, x1_x2[:, 1]), p))
    return sum(exp.(log_prob .+ log.(fisher_information_integrand.(p, x1_x2[:, 1] + x1_x2[:, 2], n.(design, x1_x2[:, 1])))))
end




mutable struct Design <: AbstractDesign
    n2::Vector{Int}
    c2::Vector{Real}
    function Design(n2::Vector{N2}, c2::Vector{C2}) where {N2<:Integer, C2<:Real}
        !all(0 .<= n2) ? error("n must be positive") : nothing
        # ToDo: check for ineffective stopping (n2 > 0 and early stopping)
        length(n2) != length(c2) ? error("n and c must be of equal length") : nothing
        new(n2, c2)
    end
end
Design(n2, c2) = Design(convert(Vector{Integer}, n2), convert(Vector{Real}, c2))
