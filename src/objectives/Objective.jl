abstract type Objective end

# reduce Base.show() to Base.string()
Base.show(io::IO, obj::Objective) = print(io, string(obj))
Base.show(io::IO, ::MIME"application/prs.juno.inline", obj::Objective) = print(io, string(obj))
Base.string(obj::Objective) = @sprintf "ToDo: implement Base.string for %s" typeof(obj)



mutable struct ExpectedScoreObjective <: Objective
    score
    orientation::Symbol
    function ExpectedScoreObjective(score, orientation)
        @assert (orientation .== (:minimise, :maximise)) |> any "orientation must be ':minimise' or ':maximise'"
        new(score, orientation)
    end
end
minimise(score::Score) = ExpectedScoreObjective(score, :minimise)
maximise(score::Score) = ExpectedScoreObjective(score, :maximise)

function Base.string(obj::ExpectedScoreObjective)
    if obj.orientation == :minimise
        return @sprintf "argmin: expected %s" string(obj.score)
    end
    if obj.orientation == :maximise
        return @sprintf "argmax: expected %s" string(obj.score)
    end
end

evaluate(obj::ExpectedScoreObjective, args...) = evaluate(obj.score, args...)
(obj::ExpectedScoreObjective)(args...) = evaluate(obj.score, args...)
integrand_x1(obj::ExpectedScoreObjective, args...; kwargs...) = integrand_x1(obj.score, args...; kwargs...)

update!(obj::ExpectedScoreObjective, x::TI, n::TI) where {TS<:Score,TI<:Integer} = update!(obj.score, x, n)


function add!(
        JuMP_model_and_indicator_variables::Tuple,
        objective::ExpectedScoreObjective,
        problem::Problem;
        partial::Tuple{TI,TI} = (0, 0)
    ) where {TI<:Integer}
    m, ind = JuMP_model_and_indicator_variables
    sense  = (objective.orientation == :minimise) ? MOI.OptimizationSense(0) : MOI.OptimizationSense(1)
    # score integrand must include a factor pmf(x1, n1, prior) but in many
    # cases it is more effective to integrate that in the score calculation
    @objective(m, sense,
        sum(
            integrand_x1(
                objective, x1, n1, n2, c2;
                partial_stage_one = partial
            ) *
            ind[(n1, x1, n2, c2)]
            for (n1, x1, n2, c2) in grid(problem)
        )
    )
end
