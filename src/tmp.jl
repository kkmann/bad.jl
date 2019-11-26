using JuMP, Ipopt, Juniper

optimizer = Juniper.Optimizer
params = Dict{Symbol,Any}()
params[:nl_solver] = with_optimizer(Ipopt.Optimizer)

m = Model(with_optimizer(optimizer, params))

@variable(m, x, start = 0.0)
@variable(m, y, start = 0.0)

f(x, y) = exp(x*y)
register(m, :f, 2, f; autodiff = true)

@NLobjective(m, Min, f(x, y))

optimize!(m)
