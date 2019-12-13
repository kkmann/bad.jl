# standard binomial
function pmf(x::TI, n::TI, p::TR; xpartial::TI = 0, npartial::TI = 0) where {TI<:Integer,TR<:Real}
    if !(0 <= xpartial <= npartial <= n)
        throw(DomainError((xpartial, npartial, n), @sprintf "0 <= xpartial=%i <= npartial=%i <= n=%i violated" x1partial n1partial n))
    end
    deltax = x - xpartial
    deltan = n - npartial
    if deltax < 0; return 0.0 end
    if deltan < deltax; return 0.0 end
    if x > n
        throw(DomainError((x, n), @sprintf "x=%i <= n=%i violated" x n))
    end
    return (gamma(deltan + 1)/gamma(deltax + 1)/gamma(deltan - deltax + 1)) *
        p^deltax * (1 - p)^(deltan - deltax)
end

function pmf(x::TI, n::TI, p::TP; xpartial::TI = 0, npartial::TI = 0) where {TI<:Integer,TP<:bad.Prior}
    expectation(
        p -> pmf(x, n, p; xpartial = xpartial, npartial = npartial),
        update(p, xpartial, npartial)
    )
end



function cdf(x::TR1, n::TI, p::TR2; xpartial::TI = 0, npartial::TI = 0) where {TI<:Integer,TR1,TR2<:Real}

    if x < 0; return 0.0 end
    if x > n; return 1.0 end
    if x == n; return 1.0 end
    if !(0 <= xpartial <= npartial <= n)
        throw(DomainError((xpartial, npartial), @sprintf "0 <= xpartial=%i <= npartial=%i <= n=%i violated" x1partial n1partial n))
    end
    deltax = x - xpartial
    deltan = n - npartial
    if deltax < 0; return 0.0 end # is x < xpartial, so is x-1 .. 0
    if deltax > deltan; return 1.0 end
    return beta_inc(deltan - deltax, deltax + 1, 1 - p, p)[1]
end

function cdf(x::TR, n::TI, p::TP; xpartial::TI = 0, npartial::TI = 0) where {TI<:Integer,TR<:Real,TP<:bad.Prior}
    expectation(
        p -> cdf(x, n, p; xpartial = xpartial, npartial = npartial),
        update(p, xpartial, npartial)
    )
end



# joint x1/x2
function pmf(x1::TI, n1::TI, x2::TI, n2::TI, p::TR; x1partial::TI = 0, n1partial::TI = 0) where {TI<:Integer,TR<:Real}
    return pmf(x2, n2, p; xpartial = 0, npartial = 0) * pmf(x1, n1, p; xpartial = x1partial, npartial = n1partial)
end

function pmf(x1::TI, n1::TI, x2::TI, n2::TI, p::TP; x1partial::TI = 0, n1partial::TI = 0) where {TI<:Integer,TP<:bad.Prior}
    expectation(
        p -> pmf(x1, n1, x2, n2, p; x1partial = x1partial, n1partial = n1partial),
        update(p, x1partial, n1partial)
    )
end

function pmf(x1::TI, x2::TI, design::TD, p::Union{TR,TP}; x1partial::TI = 0, n1partial::TI = 0) where {TI<:Integer,TR<:Real,TP<:bad.Prior,TD<:bad.AbstractDesign}
    pmf(x1, n1(design), x2, n2(design, x1), p; x1partial = x1partial, n1partial = n1partial)
end
