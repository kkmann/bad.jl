abstract type EarlyStopping end

primitive type Futility <: EarlyStopping 8 end
Futility() = reinterpret(Futility, Int8(0))

primitive type Efficacy <: EarlyStopping 8 end
Efficacy() = reinterpret(Efficacy, Int8(0))

CriticalValue = Union{Int, Futility, Efficacy}

valid(p::T) where {T<:Real} = 0 <= p <= 1
