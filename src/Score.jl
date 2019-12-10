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
(score::Score)(args...) = evaluate(score, args...)



# update the prior of a score
update!(score::TS, x::TI, n::TI) where {TS<:Score,TI<:Integer} = score.prior = update(score.prior, x, n)


function evaluate(score::TS, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR, p::TR)::TR where {TI<:Integer,TR<:Real,TS<:Score}

    error(@sprintf "ToDo: implement bad.evaluate for %s" typeof(score))
end
evaluate(score::TS, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR) where {TI<:Integer,TR<:Real,TS<:Score} = expectation(p -> evaluate(score, x1, n1, x2, n2, c2, p), score.prior)

function evaluate(score::TS, design::TD, x1::TI, x2::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    evaluate(score, x1, n1(design), x2, n2(design, x1), c2(design, x1), p)
end
evaluate(score::TS, design::TD, x1::TI, x2::TI) where {TI<:Integer,TD<:AbstractDesign,TS<:Score} = expectation(p -> evaluate(score, design, x1, x2, p), score.prior)



# conditonal expectation value given x1
function evaluate(score::TS, x1::TI, n1::TI, n2::TI, c2::TR, p::TR) where {TI<:Integer,TR<:Real,TS<:Score}

    X2 = collect(0:n2)
    return sum( evaluate.(score, x1, n1, X2, n2, c2, p) .* pmf_x2_given_x1.(X2, n2, p) )
end
evaluate(score::TS, x1::TI, n1::TI, n2::TI, c2::TR) where {TI<:Integer,TR<:Real,TS<:Score} = expectation(p -> evaluate(score, x1, n1, n2, c2, p), update(score.prior, x1, n1) )

function evaluate(score::TS, design::TD, x1::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    return evaluate(score, x1, n1(design), n2(design, x1), c2(design, x1), p)
end
evaluate(score::TS, design::TD, x1::TI) where {TI<:Integer,TD<:AbstractDesign,TS<:Score} = expectation(p -> evaluate(score, design, x1, p), update(score.prior, x1, n1(design)) )



# marginalise X1, X2
function evaluate(score::TS, design::TD, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    XX = sample_space(design)
    return sum( evaluate.(score, design, XX[:,1], XX[:,2], p) .* pdf.(XX[:,1], XX[:,2], design, p) )
end
evaluate(score::TS, design::TD) where {TI<:Integer,TD<:AbstractDesign,TS<:Score} = expectation(p -> evaluate(score, design, p), score.prior)


# compute Pr[X1=x1]*E_prior[s(x1, n1, X2, n2(x1), c2(x1))] = score_integrand(x1, n1, n2, c2)
# this is the only thing that actually used during optmisation
# feel free to implemnt specific, more efficient versions for each score!
function integrand_x1(score::TS, x1::TI, n1::TI, n2::TI, c2::TR) where {TS<:Score,TI<:Integer,TR<:Real}
    return score(x1, n1, n2, c2)*dbinom(x1, n1, score.prior)
end




mutable struct SampleSize{TP<:Prior} <: Score
    prior::TP
end
Base.string(score::SampleSize) = @sprintf "SampleSize<%s>" string(score.prior)

function evaluate(score::SampleSize, x1::TI, n1::TI, x2::TI, n2::TI, c2::TR, p::TR)::TR where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    n1 + n2
end
# no custom score_integran necessary - fallback to default





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
function integrand_x1(score::Power, x1::TI, n1::TI, n2::TI, c2::TR) where {TS<:Score,TI<:Integer,TR<:Real}
    return expectation( p-> (1 - pbinom(c2, n2, p))*dbinom(x1, n1, p), score.prior)
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
function integrand_x1(score::TypeOneErrorRate, x1::TI, n1::TI, n2::TI, c2::TR) where {TS<:Score,TI<:Integer,TR<:Real}
    return expectation( p-> (1 - pbinom(c2, n2, p))*dbinom(x1, n1, p), score.prior)
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
evaluate(score::CompositeScore, design::TD, p::TR) where {TR<:Real,TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design, p) )
evaluate(score::CompositeScore, design::TD) where {TD<:AbstractDesign} = sum( score.ω .* evaluate.(score.components, design) )

function integrand_x1(score::CompositeScore, x1::TI, n1::TI, n2::TI, c2::TR) where {TS<:Score,TI<:Integer,TR<:Real}
    return sum( score.ω .* integrand_x1.(score.components, x1, n1, n2, c2) )
end
