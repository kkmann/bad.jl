abstract type Score end
# expectation fields, must be mutable!
# - prior::Prior

# make scores iterable
Base.iterate(score::Score, state = 0) = state > 0 ? nothing : (score, state + 1)
Base.length(score::Score) = 1

# reduce Base.show() to Base.string()
Base.show(io::IO, score::Score) = print(io, string(score))
Base.show(io::IO, ::MIME"application/prs.juno.inline", score::Score) = print(io, string(score))
Base.string(score::Score) = @sprintf "ToDo: implement Base.string for %s" typeof(score)

# make scores callable
(score::Score)(args...; kwargs...) = evaluate(score, args...; kwargs...)



# update the prior of a score
update!(score::TS, x::TI, n::TI) where {TS<:Score,TI<:Integer} = score.prior = update(score.prior, x, n)



function evaluate(score::TS, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR, p::TR)::TR where {TI<:Integer,TR<:Real,TS<:Score}

    error(@sprintf "ToDo: implement bad.evaluate for %s" typeof(score))
end
function evaluate(score::TS, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR) where {TI<:Integer,TR<:Real,TS<:Score}
    expectation(
        p -> evaluate(score, x1, n1, x2, n2, c2, p),
        update(score.prior, x1 + x2, n1 + n2)
    )
end

function evaluate(score::TS, design::TD, x1::TI, x2::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    evaluate(score, x1, n1(design), x2, n2(design, x1), c2(design, x1), p)
end
function evaluate(score::TS, design::TD, x1::TI, x2::TI) where {TI<:Integer,TD<:AbstractDesign,TS<:Score}
    expectation(
        p -> evaluate(score, design, x1, x2, p),
        update(score.prior, x1 + x2, n1(design) + n2(design, x1))
    )
end



# conditonal expectation given x1
function evaluate(score::TS, x1::TI, n1::TI, n2::TI, c2::TR, p::TR) where {TI<:Integer,TR<:Real,TS<:Score}

    X2 = collect(0:n2)
    return sum( evaluate.(score, x1, n1, X2, n2, c2, p) .* pmf.(X2, n2, p) )
end
function evaluate(score::TS, x1::TI, n1::TI, n2::TI, c2::TR) where {TI<:Integer,TR<:Real,TS<:Score}
    expectation(
        p -> evaluate(score, x1, n1, n2, c2, p),
        update(score.prior, x1, n1)
    )
end

function evaluate(score::TS, design::TD, x1::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    return evaluate(score, x1, n1(design), n2(design, x1), c2(design, x1), p)
end
function evaluate(score::TS, design::TD, x1::TI) where {TI<:Integer,TD<:AbstractDesign,TS<:Score}
    expectation(
        p -> evaluate(score, design, x1, p),
        update(score.prior, x1, n1(design))
    )
end


# marginalise X1, X2
function evaluate(score::TS, design::TD, p::TR; x1partial::TI = 0, n1partial::TI = 0) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    @assert 0 <= x1partial <= n1partial <= n1(design)

    XX = sample_space(design)
    return sum( evaluate.(score, design, XX[:,1], XX[:,2], p) .* pmf.(XX[:,1], XX[:,2], design, p; x1partial = x1partial, n1partial = n1partial) )
end
function evaluate(score::TS, design::TD; x1partial::TI = 0, n1partial::TI = 0) where {TI<:Integer,TD<:AbstractDesign,TS<:Score}

    expectation(
        p -> evaluate(score, design, p; x1partial = x1partial, n1partial = n1partial),
        update(score.prior, x1partial, n1partial)
    )
end


# compute Pr[X1=x1]*E_prior[s(x1, n1, X2, n2(x1), c2(x1))] = score_integrand(x1, n1, n2, c2)
# this is the only thing that actually used during optmisation
# feel free to implement specific, more efficient versions for each score!
function integrand_x1(score::TS, x1::TI, n1::TI, n2::TI, c2::TR; x1partial::TI = 0, n1partial::TI = 0) where {TS<:Score,TI<:Integer,TR<:Real}
    score(x1, n1, n2, c2) * pmf(x1, n1, score.prior; xpartial = x1partial, npartial = n1partial)
end




mutable struct SampleSize{TP<:Prior} <: Score
    prior::TP
end
Base.string(score::SampleSize) = @sprintf "SampleSize<%s>" string(score.prior)

function evaluate(score::SampleSize, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR, p::TR)::TR where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    n1 + n2
end
# no custom score_integrand necessary - fallback to default





mutable struct Power{TP<:Prior,TR<:Real} <: Score
    prior::TP
    pmcr::TR
    Power{TP,TR}(prior::TP, pmcr::TR) where{TP<:Prior,TR<:Real} = new(pmcr <= prior, pmcr)
end
Power(prior::Prior; pmcr::Real = bounds(prior)[1]) = Power{typeof(prior),typeof(pmcr)}(prior, pmcr)
Base.string(score::Power) = @sprintf "Power<%s>" string(score.prior)

function evaluate(score::Power, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR, p::TR)::TR where {TI<:Integer,TR<:Real}
    (p >= score.pmcr) ? (x2 > c2) : 0.0
end

# since we integrate power, it is easier to do that in one wash, also
# closed form cdf is available, note that the score prior is already properly conditioned
function integrand_x1(score::Power, x1::TI, n1::TI, n2::TI, c2::TR; x1partial::TI = 0, n1partial::TI = 0) where {TS<:Score,TI<:Integer,TR<:Real}

    expectation( # conditional expected power
        p -> 1 - cdf(c2, n2, p),
        update(score.prior, x1, n1)
    ) * pmf( # conditional pmf given partial observations
        x1, n1, score.prior; xpartial = x1partial, npartial = n1partial
    )
end





mutable struct TypeOneErrorRate{TP<:Prior,TR<:Real} <: Score
    prior::TP
    pnull::TR
    TypeOneErrorRate{TP,TR}(prior::TP, pnull::TR) where{TP<:Prior,TR<:Real} = new(prior <= pnull, pnull)
end
TypeOneErrorRate(prior::Prior; pnull::Real = bounds(prior)[2]) = TypeOneErrorRate{typeof(prior),typeof(pnull)}(prior, pnull)
Base.string(score::TypeOneErrorRate) = @sprintf "TypeOneErrorRate<%s>" string(score.prior)

function evaluate(score::TypeOneErrorRate, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR, p::TR)::TR where {TI<:Integer,TR<:Real}
    (p <= score.pnull) ? (x2 > c2) : 0.0
end

# since we integrate power, it is easier to do that in one wash, also
# closed form cdf is available, note that the score prior is already properly conditioned
function integrand_x1(score::TypeOneErrorRate, x1::TI, n1::TI, n2::TI, c2::TR; x1partial::TI = 0, n1partial::TI = 0) where {TS<:Score,TI<:Integer,TR<:Real}

    expectation( # conditional expected power
        p -> 1 - cdf(c2, n2, p),
        update(score.prior, x1, n1)
    ) * pmf( # conditional pmf given partial observations
        x1, n1, score.prior; xpartial = x1partial, npartial = n1partial
    )
end




mutable struct CompositeScore <: Score
    ω
    components
end
Base.string(score::CompositeScore) = "CompositeScore"

Base.convert(CompositeScore, score::Score) = CompositeScore([1.0], [score])
Base.convert(::Type{T}, score::T)  where {T<:Score} = score



*(ω::Real, score::Score) = CompositeScore([ω], [score])
*(ω::Real, score::CompositeScore) = CompositeScore(ω .* score.ω, score.components)
function +(score1::T1, score2::T2) where {T1,T2<:Score}

    score1, score2 = convert.(CompositeScore, (score1, score2))
    return CompositeScore(vcat(score1.ω, score2.ω), vcat(score1.components, score2.components))
end
function -(score1::T1, score2::T2) where {T1,T2<:Score}

    score1, score2 = convert.(CompositeScore, (score1, score2))
    return CompositeScore(vcat(score1.ω, -score2.ω), vcat(score1.components, score2.components))
end
function update!(score::CompositeScore, x::TI, n::TI) where {TS<:Score,TI<:Integer}
    map(component -> update!(component, x, n), score.components)
end

evaluate(score::CompositeScore, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR, p::TR) where {TI<:Integer,TR<:Real} = sum( score.ω .* evaluate.(score.components, x1, n1, x2, n2, c2, p) )
evaluate(score::CompositeScore, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR) where {TI<:Integer,TR<:Real} = sum( score.ω .* evaluate.(score.components, x1, n1, x2, n2, c2) )
evaluate(score::CompositeScore, design::TD, x1::TI, x2::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design, x1, x2, p) )
evaluate(score::CompositeScore, design::TD, x1::TI, x2::TI) where {TI<:Integer,TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design, x1, x2, c2) )
evaluate(score::CompositeScore, x1::TI, n1::TI, n2::TI, c2::TR, p::TR) where {TI<:Integer,TR<:Real} = sum( score.ω .* evaluate.(score.components, x1, n1, n2, c2, p) )
evaluate(score::CompositeScore, x1::TI, n1::TI, n2::TI, c2::TR) where {TI<:Integer,TR<:Real} = sum( score.ω .* evaluate.(score.components, x1, n1, n2, c2) )
evaluate(score::CompositeScore, design::TD, x1::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design, x1, p) )
evaluate(score::CompositeScore, design::TD, x1::TI) where {TI<:Integer,TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design, x1) )
evaluate(score::CompositeScore, design::TD, p::TR; x1partial::TI = 0, n1partial::TI = 0) where {TR<:Real,TI<:Integer,TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design, p; xpartial = x1partial, npartial = n1partial) )
evaluate(score::CompositeScore, design::TD; x1partial::TI = 0, n1partial::TI = 0) where {TI<:Integer,TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design; x1partial = x1partial, n1partial = n1partial) )

function integrand_x1(score::CompositeScore, x1::TI, n1::TI, n2::TI, c2::TR; x1partial::TI = 0, n1partial::TI = 0) where {TS<:Score,TI<:Integer,TR<:Real}
    return sum( score.ω .* integrand_x1.(score.components, x1, n1, n2, c2; x1partial = x1partial, n1partial = x1partial) )
end
