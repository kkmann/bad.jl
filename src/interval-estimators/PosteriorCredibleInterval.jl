struct PosteriorCredibleInterval{TP<:Prior,TD<:AbstractDesign,TI<:Integer,TR<:Real} <: IntervalEstimator
    prior::TP
    design::TD
    XX::Array{TI,2}
    bounds::Array{TR,2}
    alpha::TR
end

string(ci::PosteriorCredibleInterval) = @sprintf "PosteriorCredibleInterval<%s>" string(ci.prior)


function PosteriorCredibleInterval(prior::Prior, design::AbstractDesign, α::Real)

    XX = sample_space(design)
    nn = size(XX, 1)

    bounds = zeros(Float64, nn, 2)
    for i in 1:nn
        x1, x2    = XX[i,:]
        posterior = update(prior, x1 + x2, n(design, x1))
        phat      = (x1 + x2) / n(design, x1)
        quantile  = prob -> Roots.find_zero(p -> cdf(p, posterior) - prob, phat)
        l = try
                quantile(α)
            catch e
                isa(e, Roots.ConvergenceFailed) ? 0.0 : e
        end
        u = try
                quantile(1 - α)
            catch e
                isa(e, Roots.ConvergenceFailed) ? 1.0 : e
        end
        u <= l ? error("credible interval collapsed") : nothing
        bounds[i,:] .= [l, u]
    end
    return PosteriorCredibleInterval{typeof(prior),typeof(design),eltype(XX),Float64}(prior, design, XX, bounds, α)
end
