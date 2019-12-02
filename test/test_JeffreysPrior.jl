design = Problem(
    ExpectedUtility(
        Beta(5, 7),
        300, -600,
        .3
    ),
    NoTypeOneErrorRateConstraint(.2, .05),
    NoPowerConstraint(.4, .2)
) |> optimise

prior = JeffreysPrior(design)

p = .025:.01:.975

pdf.(p, prior)

using Plots
Plots.plot(p, pdf.(p, prior))

posterior = update(prior, 6, 12)

mean(posterior)

Plots.plot(p, pdf.(p, update(prior, 10, 20)))
