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

function conditional(cnstr::Constraint, bounds::Tuple{Real,Real})
    cnstr = deepcopy(cnstr)
    cnstr.conditional = bounds
    return cnstr
end


mutable struct PowerConstraint <: Constraint
    score::Power
    β::Real
    conditional::Tuple{Real,Real}
end

>=(power::Power, threshold::Real) = PowerConstraint(power, 1 - threshold, (.5, .99))

function Base.string(cnstr::PowerConstraint)

    return @sprintf "%s >= %5.1f%% (given x1: %5.1f%% - %5.1f%%)" string(cnstr.score) 100*(1 - cnstr.β) 100*cnstr.conditional[1] 100*cnstr.conditional[2]
end






mutable struct TypeOneErrorRateConstraint <: Constraint
    score::Power
    α::Real
    conditional::Tuple{Real,Real}
end

<=(power::Power, threshold::Real) = TypeOneErrorRateConstraint(power, threshold, (.001, .99))

function Base.string(cnstr::TypeOneErrorRateConstraint)

    return @sprintf "%s <= %5.1f%% (given x1: %5.1f%% - %5.1f%%)" string(cnstr.score) 100*cnstr.α 100*cnstr.conditional[1] 100*cnstr.conditional[2]
end
