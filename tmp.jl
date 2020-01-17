using bad

# construct the usual prior
pnull  = .2 # response under TAU
pmcr   = pnull + .1 # minimal clinically relevant
prior1 = Beta(mean = .35, sd = .1)  # start with an subjective prior
prior2 = update(prior1, 4, 10)      # update with phase I data
prior3 = .8*prior2 + .2*Beta(1, 1)  # robustify
prior  = prior3 <= min(2*pmcr, 1.0) # restrict to plausible range

mtoer = Power(prior  | pnull)
power = Power(prior >= pmcr)
ess   = SampleSize(prior)

α, β = .05, .2
problem = Problem(
        minimise(ess),
        mtoer <= α,
        power >= 1 - β
    )
design = optimise(problem; verbosity = 0)



prior_ = (0.8*update(prior2, 2, 10) + 0.2*Beta(1, 1)) <= 0.7
x1 = 6
partial_stage_two = 2, 6
design_ = deepcopy(design)
n2_original, c2_original = n2(design, x1), c2(design, x1)
toer = design.problem.toer.score
pow  = design.problem.power.score
α_   = toer(design, x1, partial_stage_two = partial_stage_two)
β_   = 1 - pow(design, x1, partial_stage_two = partial_stage_two)
function ctoer(n2_new, c2_new)
    design_.n2[x1 + 1] = n2_new
    design_.c2[x1 + 1] = c2_new
    toer(design_, x1, partial_stage_two = partial_stage_two)
end
pow = update(pow, prior_)
function cep(n2_new, c2_new)
    design_.n2[x1 + 1] = n2_new
    design_.c2[x1 + 1] = c2_new
    pow(design_, x1, partial_stage_two = partial_stage_two)
end

function f()
    n2_new, c2_new = partial_stage_two[2], -Inf
    a, b = ctoer(n2_new, c2_new), cep(n2_new, c2_new)
    while (a > α_) | (b < 1 - β_)
        for c2_new in vcat([-Inf], 0:(n2_new - 1), [Inf])
            a, b = ctoer(n2_new, c2_new), cep(n2_new, c2_new)
            if (a <= α_) & (b >= 1 - β_)
                return n2_new, c2_new
            end
        end
        n2_new += 1
    end
end

α_, 1 - β_
n2_new, c2_new = 25, 7.0
ctoer(n2_new, c2_new), cep(n2_new, c2_new)



xx1, nn2, cc2 = adapt(design, prior_, (2, 7), verbosity = 3)

xx1, nn2, cc2 = adapt(design, prior_, (5, 7), verbosity = 3)
