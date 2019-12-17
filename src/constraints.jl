abstract type Constraint end
# expected fields (mutable)
# - score

# reduce Base.show() to Base.string()
Base.show(io::IO, cnstr::Constraint) = print(io, string(cnstr))
Base.show(io::IO, ::MIME"application/prs.juno.inline", cnstr::Constraint) = print(io, string(cnstr))

Base.string(cnstr::Constraint) = @sprintf "ToDo: implement Base.string for %s" typeof(cnstr)

integrand_x1(cnstr::Constraint, args...; kwargs...) = integrand_x1(cnstr.score, args...; kwargs...)

update!(cnstr::Constraint, x::TI, n::TI) where {TS<:Score,TI<:Integer} = update!(cnstr.score, x, n)

function update(cnstr::Constraint, prior::Prior)
    cnstr       = deepcopy(cnstr)
    cnstr.score = update(cnstr.score, prior)
    return cnstr
end


mutable struct PowerConstraint <: Constraint
    score::Power
    β::Real
    conditional::Tuple{Real,Real}
end
function subject_to(power::Power, β::Real, conditional::Tuple{Real,Real} = (.5, .99) )
    PowerConstraint(power, β, conditional)
end

function Base.string(cnstr::PowerConstraint)

    return @sprintf "%s >= %5.1f%%" string(cnstr.score) 100*(1 - cnstr.β)
end






mutable struct TypeOneErrorRateConstraint <: Constraint
    score::TypeOneErrorRate
    α::Real
    conditional::Tuple{Real,Real}
end
function subject_to(toer::TypeOneErrorRate, α::Real, conditional::Tuple{Real,Real} = (α/5, 5*α) )
    TypeOneErrorRateConstraint(toer, α, conditional)
end

function Base.string(cnstr::TypeOneErrorRateConstraint)

    return @sprintf "%s <= %5.1f%%" string(cnstr.score) 100*cnstr.α
end
