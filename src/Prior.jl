abstract type Prior end

# make priors iterable
Base.iterate(design::Prior, state = 0) = state > 0 ? nothing : (design, state + 1)
Base.length(design::Prior) = 1

Base.show(io::IO, prior::Prior) = print(io, string(prior))
Base.show(io::IO, ::MIME"application/prs.juno.inline", prior::Prior) = print(io, string(prior))

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

expected_value(f::Function, mprior::MixturePrior) = sum( mprior.ω .* expected_value.(f::Function, mprior.priors) )

mean(mprior::MixturePrior) = sum( mprior.ω .* mean.(mprior.priors) )

predictive_pmf(x, n, mprior::MixturePrior) = sum( mprior.ω .* predictive_pmf.(x, n, mprior.priors) )

function string(mprior::MixturePrior)
    n = length(mprior.priors)
    res = ""
    for i in 1:(n - 1)
        res *= (@sprintf "%.2f*" mprior.ω[i]) * string(mprior.priors[i]) * " + "
    end
    res *= (@sprintf "%.2f*" mprior.ω[n]) * string(mprior.priors[n])
    res
end
