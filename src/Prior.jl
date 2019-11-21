abstract type Prior end

# make priors iterable
Base.iterate(design::Prior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Prior) = 1



struct WeightedPrior
    ω::Real
    prior::Prior
end
*(ω::Real, prior::Prior) = WeightedPrior(ω, prior)

# make weighted priors iterable
Base.iterate(design::WeightedPrior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::WeightedPrior) = 1

condition(wprior::WeightedPrior; low::T = prior.low, high::T = prior.high) where {T<:Real} =
    wprior.ω * condition(wprior.prior; low = low, high = high)

function update(wprior::WeightedPrior, x::Int, n::Int) where {T<:Real}
    wprior.ω * update(wprior.prior, x, n)
end

integrate(wprior::WeightedPrior, values_on_pivots) = wprior.ω * integrate(wprior.prior, values_on_pivots)

mean(wprior::WeightedPrior) = wprior.ω * mean(wprior.prior)

function predictive_pmf(x, n, wprior::WeightedPrior) where {T<:Real}
    wprior.ω * integrate(wprior.prior, dbinom.(x, n, prior.pivots))
end



struct MixturePrior <: Prior
    ω::Vector{Real}
    priors::Vector{Prior}
end

+(φ::WeightedPrior, η::WeightedPrior) = sum([φ.ω, η.ω]) == 1 ? MixturePrior([φ.ω, η.ω], [φ.prior, η.prior]) : error("weights must sum to 1")

condition(mprior::MixturePrior; low::T1 = 0., high::T1 = 1.) where {T1<:Real, T2<:Real} =
    MixturePrior(mprior.ω,  condition.(mprior.priors; low = low, high = high))

update(mprior::MixturePrior, x::Int, n::Int) = MixturePrior(mprior.ω,  update.(mprior.priors, x, n))

integrate(mprior::MixturePrior, values_on_pivots) =sum( mprior.ω .* integrate.(mprior.priors, values_on_pivots) )

mean(mprior::MixturePrior) = sum( mprior.ω .* mean.(mprior.priors) )

predictive_pmf(x, n, mprior::MixturePrior) = sum( mprior.ω .* predictive_pmf.(x, n, mprior.priors) )
