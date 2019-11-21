function optimise(model::DesignIPModel; verbosity = 3, timelimit = 300)
    x1vals      = model.params["x1vals"]
    n1vals      = model.params["n1vals"]
    n2vals      = model.params["n2vals"]
    c2vals      = model.params["c2vals"]
    valid       = model.params["valid"]
    ind         = model.vars["ind"]
    n1_selected = model.vars["n1_selected"]
    optimize!(
        model.jump_model,
        with_optimizer(GLPK.Optimizer, msg_lev = verbosity, tm_lim = 1000*timelimit)
    )
    # termination_status(model.jump_model)
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
        cntr != 1 ? println([x1 cntr]) : nothing
    end
    score = objective_value(model.jump_model)
    return OptimalDesign(n2_res, c2_res, model, score)
end
