struct Prior{T<:Real}
    pdf::Function
    cdf::Function
    low::T
    high::T
    function Prior{T}(pdf::Function, cdf, low::T, high::T) where {T<:Real}
        z, prec   = quadgk(pdf, low, high)
        pdf_(p)   = low <= p <= high ? pdf(p)/z : NaN
        valid_cdf = false
        try
            valid_cdf = 0 <= f(low) <= f(high) <= 1
        catch
        end
        if valid_cdf
            cdf_ = cdf
        else
            cdf_(p) = low <= p <= high ? quadgk(pdf_, low, p)[1] : NaN
        end
        new(pdf_, cdf_, low, high,)
    end
end

Prior(pdf::Function; low::T = 0., high::T = 1., cdf = nothing) where {T<:Real} = Prior{T}(pdf, cdf, low, high)

function BetaPrior(a::T, b::T; low::T = 0., high::T = 1.) where {T<:Real}
    ϵ = sqrt(eps())
    low_, high_ = low + ϵ, high - ϵ
    z = pbeta(high_, a, b) - pbeta(low_, a, b)
    Prior(
        p -> low_ <= p <= high_ ? dbeta(p, a, b)/z : NaN;
        cdf  = p -> low <= p <= high ? pbeta(p, a, b)/z : NaN,
        low  = low_,
        high = high_
    )
end


function condition(prior::Prior{T}; low::T = prior.low, high::T = prior.high) where {T<:Real}
    low, high = max(low, prior.low), min(high, prior.high)
    z, prec   = quadgk(prior.pdf, low, high)
    pdf_(p)   = low <= p <= high ? prior.pdf(p)/z : NaN
    cdf_(p)   = low <= p <= high ? prior.cdf(p)/z : NaN
    Prior{T}(pdf_, cdf_, low, high)
end
