adapt(problem, xx1, nn1, old_design)
    m, ind, n1_selected = get_IP_model(problem)
    @constraint(m, n1_selected[nn1:end] >= 0)
    add!(jump_model, ind, problem.toer, problem, xx1, nn1, old_design)
    add!(jump_model, ind, problem.power, problem, xx1, nn1, old_design)
end
