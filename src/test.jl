

struct Problem
    objective::Objective
    type_one_error_rate_constraint::TypeOneErrorRateConstraint
    power_constraint::PowerConstraint
    params::Dict{String,Any}
end

function Problem(objective::Objective, toer::TypeOneErrorRateConstraint, power::PowerConstraint,
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
    k::Int                       = 1)
end
