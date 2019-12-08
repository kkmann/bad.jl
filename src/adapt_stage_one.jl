function pdf_k_less_n1(x1, x2, k, x_leq_k, design, prior)
    return expectation(
        p -> binomial(n2(design, x1), x2) *
            binomial(n1 - x_leq_k, x1 - x_leq_k) *
            p^(x1 + x2 - x_leq_k) *
            (1 - p)^(n(design, x1) - x1 - x2 + x_leq_k),
        prior
    )
end


function adapt_stage_one(design::TD) where {TD<:AbstractDesign}
    design = deepcopy(design)
    design.model.n1values = design.model.n1values[design.model.n1values .>= nn1]
    length(design.model.n1values) == 0 ? error("no n1values left!") : nothing
    m, ind, n1_selected = get_IP_model(design.model)
    add!(m, ind, design.model.toer, design.model, xx1, nn1, design)
    add!(m, ind, design.model.power, design.model, xx1, nn1, design)
    optimize!(m, with_optimizer(GLPK.Optimizer, msg_lev = 3, tm_lim = 1000*60))
    # termination_status(model.jump_model)
    vals = value.(ind)
    # find n1
    n1 = n1vals(design.model)[findfirst(value.(n1_selected).data .== 1.0)]
    # extract solution
    c2_res = repeat([Inf], n1 + 1)
    n2_res = zeros(n1 + 1)
    for x1 in x1vals(n1, design.model)
        cntr = 0 # make sure we have a proper solutipon
        for n2 in n2vals(n1, x1, design.model), c2 in c2vals(n1, x1, n2, design.model)
            if vals[n1, x1, n2, c2] == 1
                c2_res[x1 + 1] = c2
                n2_res[x1 + 1] = n2
                cntr          += 1
            end
        end
        cntr != 1 ? error() : nothing
    end
    score = objective_value(m)
    return OptimalDesign(n2_res, c2_res, design.model, score, design.info)
end
