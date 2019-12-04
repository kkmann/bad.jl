struct CompatibleMLE{TI<:Integer,TR<:Real,TD<:AbstractDesign} <: Estimator
    design::TD
    sample_space::Array{TI,2}
    estimates::Array{TR,1}
end

string(estimator::CompatibleMLE) = "CompatibleMLE"

function (estimator::CompatibleMLE)(x1::TI, x2::TI, design::TD) where {TI<:Integer,TD<:AbstractDesign}
    for i in 1:size(estimator.sample_space, 1)
        if [x1, x2] == estimator.sample_space[i,:]
            return estimator.estimates[i]
        end
    end
    error("(x1,x2) not found in sample space, valid observation?")
end

function CompatibleMLE(design::TD; ϵ = 1e-4, b = 1., max_iter = 10^5) where {TD<:AbstractDesign}

    XX = sample_space(design)
    # precompute test decisions and get indices for rejection/non-rejection
    # regions into XX
    decisions       = reject_null.(XX[:, 1], XX[:, 2], design)
    inds_reject     = findall(decisions .== 1)
    inds_not_reject = findall(decisions .== 0)
    # precompute mles
    mle  = MaximumLikelihoodEstimator()
    mles = mle.(XX[:, 1], XX[:, 2], design)
    nn   = size(XX, 1)
    # define versions for automatic differentiation, ll = loglikelihood
    dbinom(x, n, p)  = gamma(n + 1)/gamma(x + 1)/gamma(n - x + 1)*(1 - p)^(n - x)*p^x
    likelihood(p, i) = dbinom(XX[i,1], n1(design), p) * dbinom(XX[i,2], n2(design,XX[i,1]), p)
    smoothmax(x, b)  = sum( x .* exp.(b.*x) ./ sum(exp.(b.*x)) )
    # objective
    f(x...) = smoothmax([
        ( likelihood(x[i], i) - likelihood(mles[i], i) )^2
        for i in 1:length(x)
    ], b)
    # set up JuMP model
    m = Model(
        with_optimizer(Ipopt.Optimizer,
            max_iter        = max_iter,
            constr_viol_tol = ϵ/10
        )
    )
    # define all variables upfront, start with estimates
    @variable(m, 0 <= phat[i = 1:nn] <= 1, start = mles[i])
    # auxiliary variables for compatibility constraint
    @variable(m, a)
    @variable(m, b)
    # define compatibility constraints
    @constraint(m, a_constraint[i = inds_reject], phat[i] >= a)
    @constraint(m, b_constraint[i = inds_not_reject], phat[i] <= b)
    @constraint(m, a >= b + ϵ)
    # define stage-wise consistency constraints
    for i in 1:nn, j in 1:nn
        x1, x2   = XX[i, :]
        x1_, x2_ = XX[j, :]
        # if both are early stopping, and one has more responses we want that
        # it has a larger estimate
        if (n2(design, x1) == 0) & (n2(design, x1_) == 0) & (x1 > x1_)
                mles[i] <= mles[j] ? println([x1, x2]) : nothing
                @constraint(m, phat[i] >= phat[j] + ϵ)
        end
        # we have same x1 (=> same overall sample size),
        # estimates should be orderd by x2
        if (x1 == x1_) & (x2 > x2_)
                mles[i] <= mles[j] ? println([x1, x2]) : nothing
                @constraint(m, phat[i] >= phat[j] + ϵ)
        end
    end
    # define objective
    register(m, :f, nn, f, autodiff = true)
    # @NLobjective(m, Min, f(phat...))
    @objective(m, Min, sum( (phat .- mles).^2 ) )
    optimize!(m)
    estimates = value.(phat)

    # build object and return
    return CompatibleMLE{eltype(XX),eltype(estimates),typeof(design)}(
        design, XX, estimates
    )
end
