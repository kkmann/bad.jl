mutable struct DesignIPModel
    jump_model::Model
    vars::Dict{String,Any}
    params::Dict{String,Any}
end

function DesignIPModel(
    prior::Prior,
    p0::Real,
    α::Real,
    β::Real;
    multiple::Real               = 2,
    mrv::Real                    = p0,
    one_stage::Bool              = false,
    nmax::Int                    = one_stage ? guess_nmax(prior, mrv, p0, α, β; multiple = 2) : guess_nmax(prior, mrv, p0, α, β; multiple = multiple),
    n1min::Int                   = one_stage ? Int(round(max(.9*nmax, 5))) : Int(round(max(nmax/10, 5))),
    n1max::Int                   = one_stage ? nmax : Int(round(max(n1min, 2*nmax/3))),
    max_rel_increase::Real       = 3.,
    min_rel_increase::Real       = 1.1,
    group_sequential::Bool       = false,
    min_conditional_power::Real  = .5,
    k::Int                       = 1
)

    nmax > 200 ? error("nmax > 200, consider a continuous approximation using the R package 'adoptr'.") : nothing

    cprior       = condition(prior, low = mrv)
    n1vals       = collect(n1min:n1max)
    x1vals       = collect(0:n1max)
    n2vals       = one_stage ? [0] : collect(0:(nmax - n1max))
    c2vals_cont  = one_stage ? [] : collect(0:(nmax - n1max - 1))
    c2vals_stop  = [EarlyFutility, EarlyEfficacy]
    c2vals       = vcat(c2vals_stop, c2vals_cont) |>
        x -> convert(Vector{Union{CriticalValue, Int}}, x)

    # convenience, check whether configuration is feasible (exploit sparsity!)
    cprior = condition(prior; low = mrv)
    valid(x1, n1, n2, c2) = valid_(x1, n1, n2, c2, nmax, p0, α, cprior, n1min, n1max, max_rel_increase, min_rel_increase, mrv, min_conditional_power, k, one_stage)
    m = Model()
    @variable(m, # ind[x1, n1, n2, c2] == 1 iff n1 = n1, n2(x1) = n2, c2(x1) = c2
        ind[x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals; valid(x1, n1, n2, c2)],
        Bin
    )
    # need to make sure that exactly one n1 value is selected, IP or trick via
    # auxiliary variables: n1_selected[n1] == 1 iff n1 = n1
    @variable(m, n1_selected[n1 in n1vals], Bin)
    @constraint(m, sum(n1_selected[n1] for n1 in n1vals) == 1)
    if group_sequential # same for n2 in group sequential case
        @variable(m, n2_selected[n2 in n2vals], Bin)
        @constraint(m, sum(n2_selected[n2] for n2 in n2vals) == 1)
    end
    # implement all constraints conditional on n1
    for n1 in n1vals
        # make sure that n1_selected[n1] is 1 iff any of the other values is assigned
        @constraint(m,
            2*nmax * n1_selected[n1] >= sum(ind[x1, n1, n2, c2] for x1 in x1vals, n2 in n2vals, c2 in c2vals if valid(x1, n1, n2, c2))
        )
        # make sure that n2_selected[n2] is 1 iff any of the other values is assigned
        if group_sequential
            for n2 in n2vals
                @constraint(m,
                    2*nmax * n2_selected[n2] >= sum(ind[x1, n1, n2, c2] for x1 in x1vals, c2 in c2vals_cont if valid(x1, n1, n2, c2))
                )
            end
        end
        for x1 in 0:n1
            # make its a function in x1
            @constraint(m,
                n1_selected[n1] == sum(ind[x1, n1, n2, c2] for n2 in n2vals, c2 in c2vals if valid(x1, n1, n2, c2))
            )
            # make sure that n1(x1) is constant
            if x1 > 0
                @constraint(m,
                    sum(ind[x1,     n1, n2, c2] for n2 in n2vals for c2 in c2vals if valid(x1, n1, n2, c2) & valid(x1 - 1, n1, n2, c2)) ==
                    sum(ind[x1 - 1, n1, n2, c2] for n2 in n2vals for c2 in c2vals if valid(x1, n1, n2, c2) & valid(x1 - 1, n1, n2, c2))
                )
            end
            # make sure we have contiguous stopping for efficacy / futility
            if ((x1 > 0) & all(valid.([x1 - 1, x1], n1, 0, [EarlyFutility])))
                @constraint(m,
                    ind[x1 - 1, n1, 0, EarlyFutility] >= ind[x1, n1, 0, EarlyFutility]
                )
            end
            if ((x1 < maximum(n1vals)) & all(valid.([x1, x1 + 1], n1, 0, [EarlyEfficacy])))
                @constraint(m,
                    ind[x1 + 1, n1, 0, EarlyEfficacy] >= ind[x1, n1, 0, EarlyEfficacy]
                )
            end
        end
    end
    # maximal type one error rate constraint, TODO: externalize?
    @constraint(m,
        sum(
            dbinom(x1, n1, p0) * power(n2, c2, p0) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        ) <= α
    )
    params = Dict(
        "prior"                 => prior,
        "p0"                    => p0,
        "α"                     => α,
        "β"                     => β,
        "mrv"                   => mrv,
        "min_conditional_power" => min_conditional_power,
        "one_stage"             => one_stage,
        "group_sequential"      => group_sequential,
        "nmax"                  => nmax,
        "n1min"                 => n1min,
        "n1max"                 => n1max,
        "max_rel_increase"      => max_rel_increase,
        "min_rel_increase"      => min_rel_increase,
        "multiple"              => multiple,
        "k"                     => k,
        "cprior"                => cprior,
        "n1vals"                => n1vals,
        "x1vals"                => x1vals,
        "n2vals"                => n2vals,
        "c2vals_cont"           => c2vals_cont,
        "c2vals_stop"           => c2vals_stop,
        "c2vals"                => c2vals,
        "cprior"                => cprior,
        "valid"                 => valid
    )
    vars = Dict(
        "ind"         => ind,
        "n1_selected" => n1_selected
    )
    return DesignIPModel(m, vars, params)

end

function valid_(x1, n1, n2, c2, nmax, p0, α, cprior, n1min, n1max, max_rel_increase, min_rel_increase, mrv, min_conditional_power, k, one_stage)

    x1 < 0 ?
        (return false) : nothing
    x1 > n1 ?
        (return false) : nothing
    one_stage & (n2 > 0) ?
        (return false) : nothing
    n1 + n2 > nmax ?
        (return false) : nothing
    (n1 + n2)/n1 > max_rel_increase ?
        (return false) : nothing
    # minimum relative increae only relevant when not stopping!
    ( !early_stop(c2) & ((n1 + n2)/n1 < min_rel_increase) ) ?
        (return false) : nothing
    # c2 can only be larger or equal to n2 if it is infinite = EarlyFutility
    ( (c2 >= n2) & (c2 != EarlyFutility) ) ?
        (return false) : nothing
    # ensure n2 is zero when stopping early
    ( !early_stop(c2) & (n2 == 0) ) ?
        (return false) : nothing
    ( early_stop(c2) & (n2 > 0) ) ?
        (return false) : nothing
    # curtailment: if we can reject based on stage one data alone, do so!
    # must use k >=1 to account for error inflation by using two stage design!
    ( (x1 > findfirst(1 .- pbinom.(0:n1, n1, p0) .<= α) + k) & (n2 > 0) ) ?
        (return false) : nothing
    # check conditional power constraint
    if min_conditional_power > 0 # save time if not!
        ( !early_stop(c2) & (power(x1, n1, n2, c2, cprior) <= min_conditional_power) ) ?
                (return false) : nothing
    end
    return true # noting to object, return true
end
