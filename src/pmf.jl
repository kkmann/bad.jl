# standard binomial
function pmf(x::TI, n::TI, p::TR; partial::Tuple{TI,TI} = (0, 0)) where {TI<:Integer,TR<:Real}

    if !(0 <= partial[1] <= partial[2] <= n)
        throw(DomainError((partial..., n), @sprintf "0 <= xpartial=%i <= npartial=%i <= n=%i violated" partial[1] partial[2] n))
    end
    deltax = x - partial[1]
    deltan = n - partial[2]
    if deltax < 0; return 0.0 end
    if deltan < deltax; return 0.0 end
    if x > n
        throw(DomainError((x, n), @sprintf "x=%i <= n=%i violated" x n))
    end
    return (gamma(deltan + 1)/gamma(deltax + 1)/gamma(deltan - deltax + 1)) *
        p^deltax * (1 - p)^(deltan - deltax)
end

function pmf(x::TI, n::TI, p::TP; partial::Tuple{TI,TI} = (0, 0)) where {TI<:Integer,TP<:Prior}
    expectation(
        p -> pmf(x, n, p; partial = partial),
        update(p, partial[1], partial[2])
    )
end



function cdf(x::TR1, n::TI, p::TR2; partial::Tuple{TI,TI} = (0, 0)) where {TI<:Integer,TR1,TR2<:Real}

    if x < 0; return 0.0 end
    if x > n; return 1.0 end
    if x == n; return 1.0 end
    if !(0 <= partial[1] <= partial[2] <= n)
        throw(DomainError((partial..., n), @sprintf "0 <= xpartial=%i <= npartial=%i <= n=%i violated" partial[1] partial[2] n))
    end
    deltax = x - partial[1]
    deltan = n - partial[2]
    if deltax < 0; return 0.0 end # is x < xpartial, so is x-1 .. 0
    if deltax > deltan; return 1.0 end
    return beta_inc(deltan - deltax, deltax + 1, 1 - p, p)[1]
end

function cdf(x::TR, n::TI, p::TP; partial::Tuple{TI,TI} = (0, 0)) where {TI<:Integer,TR<:Real,TP<:bad.Prior}
    expectation(
        p -> cdf(x, n, p; partial = partial),
        update(p, partial[1], partial[2])
    )
end
