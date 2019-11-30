problem = Problem(
    ExpectedUtility(
        Beta(5, 7),
        300, -150,
        .3
    ),
    NoTypeOneErrorRateConstraint(.2, .05),
    NoPowerConstraint(.4, .2)
)

design = optimise(problem)

plot(design)

import Gadfly
p = 0:.01:1
Gadfly.plot(x = collect(p), y = power.(design, p))

power(design, .2)
power(design, .3)
power(design, .4)
