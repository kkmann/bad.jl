using Test

p0    = .2
pmcr  = .3
prior = Beta(5, 7)
α, β  = .05, .2

design = Problem(
        minimise_expected_sample_size(prior),
        maximal_type_one_error_rate(p0, α),
        minimal_expected_power(prior, pmcr, 1 - β),
    ) |> optimise

samplesize = SampleSize(prior)

XX = sample_space(design)

pdf.(XX[:,1], XX[:,2], design, prior) |> sum

evaluate.(samplesize, design, XX[:,1], XX[:,2], .5)


evaluate(samplesize, design, 0, 0, .2)
evaluate(samplesize, design, 0, 0)
evaluate(samplesize, design, 0, .2)
evaluate(samplesize, design, 0)
evaluate(samplesize, design, .2)
evaluate(samplesize, design)

@time evaluate(samplesize, design)

@time expectation(p -> expected_sample_size(design, p), prior)



pow = Power(prior, pmcr)

evaluate(pow, design, 0, 0, .2)
evaluate(pow, design, 0, 0)
evaluate(pow, design, 6, 0)
evaluate(pow, design, 6, .2)
evaluate(pow, design, 5)
evaluate(pow, design, .2)
evaluate(pow, design)

toer = TypeOneErrorRate(prior, p0)
evaluate(toer, design, 0, 0, .2)
evaluate(toer, design, 0, 0)
evaluate(toer, design, 6, 0)
evaluate(toer, design, 6, .2)
evaluate(toer, design, 5)
evaluate(toer, design, .2)
evaluate(toer, design)
