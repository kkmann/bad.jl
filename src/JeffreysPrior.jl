
# """
#     jeffreysprior(design::Design)
# Computes the Jeffreys prior of any given design.
# # Parameters
# | Parameter    | Description |
# | -----------: | :---------- |
# | design       | a Design |
# # Return Value
# A function with signature `prior{T<:Real}(p::T)::Real` where `p` is the response
# probability and `prior(p)` the PDF of the Jeffres prior at `p`.
# # Examples
# ```julia-repl
# julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
# julia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)
# julia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())
# julia> f = jeffreysprior(design)
# julia> f(.5)
# ```
# """
# function jeffreysprior(design::Design)
#
#   function sqrtfi(p::Float64)::Float64
#
#     res  = 0.0
#     n1   = interimsamplesize(design)
#     nmax = maximum(samplesize(design))
#     supp = support(design)
#     for i in 1:size(supp, 1)
#         x1, x2 = supp[i, 1], supp[i, 2]
#         x   = x1 + x2
#         n   = samplesize(design, x1)
#         res += binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2))*BigFloat(p^x*(1 - p)^(n - x)*(x/p - (n - x)/(1 - p))^2)
#     end
#     return sqrt(res)
#
#   end
#
#   z = QuadGK.quadgk(sqrtfi, 0, 1, atol    = 0.001)[1] # exact integration from 0 to 1 is expensive!
#
#   function prior{T<:Real}(p::T)::Real
#
#     @checkprob p
#     return sqrtfi(p)/z
#
#   end
#
#   return prior
#
# end
