using Test

design = Problem(
        ExpectedUtility(
            Beta(5, 7),
            300, -600,
            .3
        ),
        NoTypeOneErrorRateConstraint(.2, .05),
        NoPowerConstraint(.4, .2)
    ) |>
    pr -> optimise(pr, mle_incompatible_is_error = false)

prior = JeffreysPrior(design)

p = .025:.01:.975

pdf.(p, prior)

posterior = update(prior, 6, 12)

mean(posterior)
