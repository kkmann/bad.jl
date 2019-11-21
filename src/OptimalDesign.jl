mutable struct OptimalDesign <: AbstractDesign
    n2::Vector{Int}
    c2::Vector{CriticalValue}
    model::DesignIPModel
    score::Real
    # check n2/c2
end
