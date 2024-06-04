abstract type Basis end

struct LaurentPolynomialBasis{T} <: Basis
    elements::Vector{LaurentPolynomial{T}}
end

@inline Base.size(lpb::LaurentPolynomialBasis) = size(lpb.elements)
@inline Base.getindex(lpb::LaurentPolynomialBasis, n::Int) =  lpb.elements[n] 
@inline Base.eachindex(lpb::LaurentPolynomialBasis) = eachindex(lpb.elements)


# Mass matrix
function mass_matrix(lpb::LaurentPolynomialBasis)
    A = [scalar_product(lpb[i], lpb[j]) for i in eachindex(lpb) for j in eachindex(lpb)]
    reshape(A, (size(lpb)[1],size(lpb)[1]))
end

function mass_matrix(lpb::LaurentPolynomialBasis, a::Real, b::Real)
    A = [scalar_product(lpb[i], lpb[j], a, b) for i in eachindex(lpb) for j in eachindex(lpb)]
    reshape(A, (size(lpb)[1],size(lpb)[1]))
end

function weight_mass_matrix(lpb::LaurentPolynomialBasis, weight::LaurentPolynomial)
    A = [scalar_product(weight * lpb[i], lpb[j]) for i in eachindex(lpb) for j in eachindex(lpb)]
    reshape(A, (size(lpb)[1],size(lpb)[1]))
end

function weight_mass_matrix(lpb::LaurentPolynomialBasis, weight::LaurentPolynomial, a::Real, b::Real)
    A = [scalar_product(weight * lpb[i], lpb[j], a, b) for i in eachindex(lpb) for j in eachindex(lpb)]
    reshape(A, (size(lpb)[1],size(lpb)[1]))
end

function weight_mass_matrix(lpb::LaurentPolynomialBasis, n::Int, a::Real, b::Real)
    weight_mass_matrix(lpb::LaurentPolynomialBasis, Monomial(n), a::Real, b::Real)
end