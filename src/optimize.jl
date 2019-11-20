function get_optimal_design(
    nmax::Int,
    p0::Real,
    α::Real,
    β::Real,
    prior::Prior;
    n1min::Int                   = Int(round(max(nmax/10, 5))),
    n1max::Int                   = Int(round(max(n1min, 2*nmax/3))),
    max_rel_increase::Real       = 3.,
    min_rel_increase::Real       = 1.1,
    mrv::Real                    = p0,
    min_conditional_power::Real  = 0,
    k::Int                       = 1,
    max_seconds::Int             = 60,
    verbose::Int                 = 3
)

    p1 = .4

    n1vals       = collect(n1min:n1max)
    x1vals       = collect(0:n1max)
    n2vals       = collect(0:(nmax - n1max))
    c2vals_cont  = collect(0:(nmax - n1max - 1))
    c2vals_stop  = [EarlyFutility, EarlyEfficacy]
    c2vals       = vcat([EarlyFutility], 0:(maximum(n2vals) - 1), [EarlyEfficacy]) |>
        x -> convert(Vector{Union{CriticalValue, Int}}, x)
    cprior       = condition(prior, low = mrv)
    # convenience, check whether configuration is feasible (exploit sparsity!)
    valid(x1, n1, n2, c2) = valid_(x1, n1, n2, c2, nmax, p0, α, prior, n1min, n1max, max_rel_increase, min_rel_increase, mrv, min_conditional_power, k)
    # define model basis with indicator variables for all possible feasible configurations
    m = Model()
    println("creating indicator variables...")
    @variable(m, # ind[x1, n1, n2, c2] == 1 iff n1 = n1, n2(x1) = n2, c2(x1) = c2
        ind[x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals; valid(x1, n1, n2, c2)],
        Bin
    )
    @printf "IP model with %i variables created.\n\r" length(ind)
    # need to make sure that exactly one n1 value is selected, IP or trick via
    # auxiliary variables: n1_selected[n1] == 1 iff n1 = n1
    @variable(m, n1_selected[n1 in n1vals], Bin)
    @constraint(m, sum(n1_selected[n1] for n1 in n1vals) == 1)
    println("defining feasibility constraints...")
    # implement all constraints conditional on n1
    for n1 in n1vals
        # make sure that n1_selected[n1] is 1 iff any of the other values is assigned
        @constraint(m,
            nmax * n1_selected[n1] >= sum(ind[x1, n1, n2, c2] for x1 in x1vals, n2 in n2vals, c2 in c2vals if valid(x1, n1, n2, c2))
        )
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
    println("done.")
    println("defining type one error rate constraint...")
    # maximal type one error rate constraint
    @constraint(m,
        sum(
            dbinom(x1, n1, p0) * power(x1, n2, c2, p0) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        ) <= α
    )
    println("done.")
    println("defining power constraint...")
    # power constraint
    @constraint(m,
        sum(
            dbinom(x1, n1, p1) * power(x1, n2, c2, p1) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        ) >= 1 - β
    )
    println("done.")
    println("defining objective...")
    # objective: minimize expected sample size
    @objective(m, Min,
        sum(
            (n1 + n2) * dbinom(x1, n1, p1) * ind[x1, n1, n2, c2] for
            x1 in x1vals, n1 in n1vals, n2 in n2vals, c2 in c2vals if
            valid(x1, n1, n2, c2)
        )
    )
    println("done.")
    optimize!(m, with_optimizer(GLPK.Optimizer, msg_lev = verbose, tm_lim = 1000*max_seconds))
    termination_status(m)
    # todo: check solution!
    objective_value(m)
    vals = value.(ind)
    # find n1
    n1 = n1vals[findfirst(value.(n1_selected).data .== 1.0)]
    # extract solution
    c2_res = convert(Vector{CriticalValue}, repeat([EarlyFutility], n1 + 1))
    n2_res = zeros(n1 + 1)
    for x1 in 0:n1
        cntr = 0 # make sure we have a proper solutipon
        for n2 in n2vals, c2 in c2vals
            if valid(x1, n1, n2, c2)
                if vals[x1, n1, n2, c2] == 1
                    c2_res[x1 + 1] = c2
                    n2_res[x1 + 1] = n2
                    cntr          += 1
                end
            end
        end
        cntr  > 1 ? println([x1 cntr]) : nothing
        cntr == 0 ? println([x1 cntr]) : nothing
    end
    # check solution
    return Design(n2_res, c2_res)
end



function valid_(x1, n1, n2, c2, nmax, p0, α, prior, n1min, n1max, max_rel_increase, min_rel_increase, mrv, min_conditional_power, k)
    x1 > n1 ?
        (return false) : nothing
    !(n1 + n2 <= nmax) ?
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
    # if min_conditional_power > 0 # save time if not!
    #     ( early_stop(c2) & (power(x1, n1, n2, c2, cprior) <= min_conditional_power) ) ?
    #         (return false) : nothing
    # end
    return true # noting to object, return true
end
