using Test, bad; import Plots



@info "testing utility functions..."
@time include("test_util.jl")

@info "testing probability mass function..."
@time include("test_pmf.jl")

@info "testing priors..."
@time include("test_Beta.jl")

@info "testing scores..."
@time include("test_scores.jl")

@info "testing point estimatos..."
@time include("test_point-estimators.jl")

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

@info "... done - stay vigilent!"
