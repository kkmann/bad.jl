using Test, bad

include("test_util.jl")

include("test_Beta.jl")
include("test_MixturePrior.jl")
include("test_JeffreysPrior.jl")

include("test_Design.jl")
include("test_MiniMaxSampleSize.jl")

include("test_MaximumLikelihoodEstimator.jl")
include("test_PosteriorMean.jl")
include("test_RaoBlackwellEstimator.jl")
include("test_CompatibleMLE.jl")

include("test_interval.jl")
