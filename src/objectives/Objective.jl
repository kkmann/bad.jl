abstract type Objective end


(::Objective)(design::AbstractDesign) = error("not implemented")
(::Objective)(design::AbstractDesign, x1::Integer) = error("not implemented")
(::Objective)(design::AbstractDesign, x1::Integer, x2::Integer) = error("not implemented")
(::Objective)(design::AbstractDesign, x1::Integer, x2::Integer, p::Real) = error("not implemented")

update!(objective::Objective, prior::Prior) = objective

# add!(jump_model, ind, objective::Objective, problem::Problem)
