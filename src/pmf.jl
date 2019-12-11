function pmf(x1::TI, n1::TI, x2::TI, n2::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    return dbinom(x1, n1, p) * dbinom(x2, n2, p)
end
function pmf(x1::TI, n1::TI, x2::TI, n2::TI, prior::TP) where {TI<:Integer,TP<:Prior,TD<:AbstractDesign}
    return expectation(p -> pmf(x1, n1, x2, n2, p), prior)
end
function pmf(design::TD, x1::TI, x2::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    return pmf(x1, n1(design), x2, n2(design, x1), p)
end
function pmf(design::TD, x1::TI, x2::TI, prior::TP) where {TI<:Integer,TP<:Prior,TD<:AbstractDesign}
    return expectation(p -> pmf(design, x1, x2, p), prior)
end

# marginal
function pmf(x1::TI, n1::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    return dbinom(x1, n1, p)
end
function pmf(x1::TI, n1::TI, prior::TP) where {TI<:Integer,TP<:Prior,TD<:AbstractDesign}
    return expectation(p -> pmf(x1, n1, p), prior)
end
function pmf(design::TD, x1::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    return pmf(x1, n1(design), p)
end
function pmf(design::TD, x1::TI, prior::TP) where {TI<:Integer,TP<:Prior,TD<:AbstractDesign}
    return expectation(p -> pmf(design, x1, p), prior)
end

# conditional
function pmf_x2_given_x1(x2::TI, n2::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    return dbinom(x2, n2, p)
end
function pmf_x2_given_x1(x2::TI, n2::TI, prior::TP) where {TI<:Integer,TP<:Prior,TD<:AbstractDesign}
    return expectation(p -> pmf_x2_given_x1(x2, n2, p), prior)
end
function pmf_x2_given_x1(design::TD, x1::TI, x2::TI, p::TR) where {TI<:Integer,TR<:Real,TD<:AbstractDesign}
    return pmf_x2_given_x1(x2, n2(design, x1), p)
end
function pmf_x2_given_x1(design::TD, x1::TI, x2::TI, prior::TP) where {TI<:Integer,TP<:Prior,TD<:AbstractDesign}
    return expectation(p -> dbinom(x2, n2(design,x1), p), prior)
end




# conditional on x(0, n1 ∧ τ], x(n1 ∧ τ, τ]
