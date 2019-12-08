abstract type Score end
# expected fields:
# - prior::Prior

# make scores iterable
Base.iterate(score::Score, state = 0) = state > 0 ? nothing : (score, state + 1)
Base.length(score::Score) = 1

# reduce Base.show() to Base.string()
Base.show(io::IO, score::Score) = print(io, string(score))
Base.show(io::IO, ::MIME"application/prs.juno.inline", score::Score) = print(io, string(score))
Base.string(score::Score) = @sprintf "ToDo: implement Base.string for %s" typeof(score)



function evaluate(score::TS, design::TD, x1::TI, x2::TI, p::TR)::TR where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    error(@sprintf "ToDo: implement bad.evaluate for %s" typeof(score))
end
# integrate over prior
function evaluate(score::TS, design::TD, x1::TI, x2::TI) where {TI<:Integer,TD<:AbstractDesign,TS<:Score}

    return expectation(p -> evaluate(score, design, x1, x2, p), score.prior)
end

# conditonal expected value given x1
function evaluate(score::TS, design::TD, x1::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    X2 = sample_space(design, x1)
    return sum( evaluate.(score, design, x1, X2, p) .* pmf_x2_given_x1.(design, x1, X2, p) )
end
function evaluate(score::TS, design::TD, x1::TI) where {TI<:Integer,TD<:AbstractDesign,TS<:Score}

    return sum( expectation(p -> evaluate(score, design, x1, p), update(score.prior, x1, n1(design)) ) )
end

# marginalise X1, X2
function evaluate(score::TS, design::TD, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign,TS<:Score}

    XX = sample_space(design)
    return sum( evaluate.(score, design, XX[:,1], XX[:,2], p) .* pdf.(XX[:,1], XX[:,2], design, p) )
end
function evaluate(score::TS, design::TD) where {TI<:Integer,TD<:AbstractDesign,TS<:Score}

    return expectation(p -> evaluate(score, design, p), score.prior)
end



struct SampleSize{TP<:Prior} <: Score
    prior::TP
end
Base.string(score::SampleSize) = @sprintf "SampleSize<%s>" string(score.prior)

function evaluate(score::SampleSize, design::TD, x1::TI, x2::TI, p::TR)::TR where {TI<:Integer,TR<:Real,TD<:AbstractDesign}

    n(design, x1)
end


struct Power{TP<:Prior,TR<:Real} <: Score
    prior::TP
    pmcr::TR
    Power{TP,TR}(prior::TP, pmcr::TR) where{TP<:Prior,TR<:Real} = new(pmcr <= prior, pmcr)
end
Power(prior::Prior, pmcr::Real) = Power{typeof(prior),typeof(pmcr)}(prior, pmcr)
Base.string(score::Power) = @sprintf "Power<%s>" string(score.prior)


function evaluate(score::Power, design::TD, x1::TI, x2::TI, p::TR)::TR where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    p >= score.pmcr ? (x2 > c2(design, x1)) : 0.0
end


struct TypeOneErrorRate{TP<:Prior,TR<:Real} <: Score
    prior::TP
    pnull::TR
    TypeOneErrorRate{TP,TR}(prior::TP, pnull::TR) where{TP<:Prior,TR<:Real} = new(prior <= pnull, pnull)
end
TypeOneErrorRate(prior::Prior, pnull::Real) = TypeOneErrorRate{typeof(prior),typeof(pnull)}(prior, pnull)
Base.string(score::TypeOneErrorRate) = @sprintf "TypeOneErrorRate<%s>" string(score.prior)


function evaluate(score::TypeOneErrorRate, design::TD, x1::TI, x2::TI, p::TR)::TR where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    p <= score.pnull ? (x2 > c2(design, x1)) : 0.0
end
