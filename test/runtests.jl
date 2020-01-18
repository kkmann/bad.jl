using Test, bad



@info "testing probability mass function..."
@time include("test_pmf.jl")

@info "testing priors..."
@time include("test_Beta.jl")

@info "testing scores..."
@time include("test_scores.jl")

@info "testing point estimators..."
@time include("test_point-estimators.jl")

@info "testing p-values..."
@time include("test_p-values.jl")

@info "testing interval estimators..."
@time include("test_intervals.jl")

@info "testing one stage designs..."
@time include("test_one_stage.jl")

@info "testing vs simons designs..."
@time include("test_simons.jl")

@info "testing minimmax sample size..."
@time include("test_MiniMaxSampleSize.jl")

@info "testing expected utility..."
@time include("test_ExpectedUtility.jl")

@info "testing design adaptations..."
@time include("test_adapt.jl")

@info "... done - stay vigilent!"
