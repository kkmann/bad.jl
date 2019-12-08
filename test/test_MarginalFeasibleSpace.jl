using Test

prior = Beta(mean = .4, sd = .1)
pnull = .2
pmcr  = .3
α, β  = .05, .2

mfs = bad.MarginalFeasibleSpace(pnull, pmcr)
