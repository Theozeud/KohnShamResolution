abstract type Basis end

struct LaurentPolynomialBasis{TL <: AbstractLaurentPolynomial} <: Basis
    elements::Vector{TL}
end

@inline Base.isempty(lpb::LaurentPolynomialBasis) = isempty(lpb.elements)
@inline Base.size(lpb::LaurentPolynomialBasis) = size(lpb.elements)
@inline Base.length(lpb::LaurentPolynomialBasis) = length(lpb.elements)
@inline Base.getindex(lpb::LaurentPolynomialBasis, n::Int) =  lpb.elements[n] 
@inline Base.eachindex(lpb::LaurentPolynomialBasis) = eachindex(lpb.elements)
@inline Base.iterate(lpb::LaurentPolynomialBasis, state = 1) = state > length(lpb) ? nothing : (lpb[state],state+1)
@inline Base.first(lpb::LaurentPolynomialBasis) = first(lpb.elements)


# Mass matrix
function mass_matrix(lpb::LaurentPolynomialBasis{TL}) where TL
    A = [scalar_product(lpb[i], lpb[j]) for i in eachindex(lpb) for j in eachindex(lpb)]
    reshape(A, (size(lpb)[1],size(lpb)[1]))
end

function mass_matrix(lpb::LaurentPolynomialBasis, a::Real, b::Real)
    T = eltype(first(lpb))
    A = zeros(T, (length(lpb), length(lpb)))
    for i ∈ eachindex(lpb)
        for j ∈ eachindex(lpb)
            if i > j
                A[i,j] = A[j,i] 
            else
                A[i,j] = scalar_product(lpb[i], lpb[j], a, b)
            end
        end
    end
    A
end

function weight_mass_matrix(lpb::LaurentPolynomialBasis, weight::LaurentPolynomial)
    A = [scalar_product(weight * lpb[i], lpb[j]) for i in eachindex(lpb) for j in eachindex(lpb)]
    reshape(A, (size(lpb)[1],size(lpb)[1]))
end

function weight_mass_matrix(lpb::LaurentPolynomialBasis, weight::LaurentPolynomial, a::Real, b::Real)
    T = eltype(first(lpb))
    A = zeros(T, (length(lpb), length(lpb)))
    for i ∈ eachindex(lpb)
        for j ∈ eachindex(lpb)
            if i > j
                A[i,j] = A[j,i] 
            else
                A[i,j] = scalar_product(weight * lpb[i], lpb[j], a, b)
            end
        end
    end
    A
end

function weight_mass_matrix(lpb::LaurentPolynomialBasis, n::Int, a::Real, b::Real)
    weight_mass_matrix(lpb::LaurentPolynomialBasis, Monomial(n), a::Real, b::Real)
end

# Deriv
function deriv!(lpb::LaurentPolynomialBasis)
    for i in eachindex(lpb)
        deriv!(lpb[i])
    end
    lpb
end

function deriv(lpb::LaurentPolynomialBasis)
    TL = eltype(lpb.elements)
    deriv_laurent = TL[]
    for i in eachindex(lpb)
        push!(deriv_laurent, deriv(lpb[i]))
    end
    LaurentPolynomialBasis(deriv_laurent)
end

# Build LaurentPolynomial on basis
function build_on_basis(basis::LaurentPolynomialBasis{TL}, coeff) where TL
    if isempty(basis)
        @error "The basis is empty."
    end
    @assert length(coeff) == length(basis)
    poly = zero(first(basis.elements))
    for i ∈ eachindex(basis)
        poly += coeff[i] * basis[i]
    end
    poly
end
