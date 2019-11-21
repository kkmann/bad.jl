abstract type Prior end

# make priors iterable
Base.iterate(design::Prior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Prior) = 1

struct MixturePrior <: Prior
    ω::Vector{Real}
    priors::Vector{Prior}
end
*(ω::Real, prior::Prior) = MixturePrior([ω], [prior])
+(φ::MixturePrior, η::MixturePrior) = MixturePrior(vcat(φ.ω, η.ω), vcat(φ.priors, η.priors))

is_proper(mprior::MixturePrior) = sum(mprior.ω) == 1

condition(mprior::MixturePrior; low::T1 = 0., high::T1 = 1.) where {T1<:Real, T2<:Real} =
    MixturePrior(mprior.ω,  condition.(mprior.priors; low = low, high = high))

update(mprior::MixturePrior, x::Int, n::Int) = MixturePrior(mprior.ω,  update.(mprior.priors, x, n))

integrate(mprior::MixturePrior, values_on_pivots) =sum( mprior.ω .* integrate.(mprior.priors, values_on_pivots) )

mean(mprior::MixturePrior) = sum( mprior.ω .* mean.(mprior.priors) )

predictive_pmf(x, n, mprior::MixturePrior) = sum( mprior.ω .* predictive_pmf.(x, n, mprior.priors) )
