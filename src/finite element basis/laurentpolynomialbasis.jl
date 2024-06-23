abstract type Basis end

abstract type AbstractLaurentPolynomialBasis end

struct LaurentPolynomialBasis{TL <: AbstractLaurentPolynomial} <: AbstractLaurentPolynomialBasis
    mesh::OneDMesh
    elements::Vector{TL}
    cross_index
    function LaurentPolynomialBasis(elements)
        mesh = first(elements).mesh
        cross_index = CartesianIndex[]
        for i in eachindex(elements)
            for j in 1:i
                if !isempty(intersect(elements[i].index, elements[j].index))
                    push!(cross_index, CartesianIndex(i, j))
                end
            end
        end
        new{eltype(elements)}(mesh, elements, cross_index)
    end
end

@inline Base.isempty(lpb::LaurentPolynomialBasis) = isempty(lpb.elements)
@inline Base.size(lpb::LaurentPolynomialBasis) = size(lpb.elements)
@inline Base.length(lpb::LaurentPolynomialBasis) = length(lpb.elements)
@inline Base.getindex(lpb::LaurentPolynomialBasis, n::Int) =  lpb.elements[n] 
@inline Base.eachindex(lpb::LaurentPolynomialBasis) = eachindex(lpb.elements)
@inline Base.iterate(lpb::LaurentPolynomialBasis, state = 1) = state > length(lpb) ? nothing : (lpb[state],state+1)
@inline Base.first(lpb::LaurentPolynomialBasis) = first(lpb.elements)
@inline bottom_type(lpb::LaurentPolynomialBasis) = eltype(first(lpb))

"""
    mass_matrix(lpb::LaurentPolynomialBasis)

Compute the mass matrix of the LaurentPolynomialBasis.

To do that, it computes the scalar product between each elements of the basis according to
lpb.cross_index that knows which index of the matrix this scalar_product is a priori not null. 
"""
function mass_matrix(lpb::LaurentPolynomialBasis)
    T = eltype(first(lpb))
    A = zeros(T, (length(lpb), length(lpb)))
    for I ∈ lpb.cross_index
        @inbounds A[I[1],I[2]] = scalar_product(lpb[I[1]], lpb[I[2]], first(lpb.mesh), last(lpb.mesh))
        @inbounds A[I[2],I[1]] = A[I[1],I[2]] 
    end
    A
end

"""
    weight_mass_matrix(lpb::LaurentPolynomialBasis, weight::LaurentPolynomial)

Compute a weighted mass matrix of the LaurentPolynomialBasis, i.e the matrix made of 
scalar product of the product of weight and an element of the basis with an element of the
basis. 

The computations are done according to lpb.cross_index that knows which index of the matrix 
this scalar_product is a priori not null. 
"""
function weight_mass_matrix(lpb::LaurentPolynomialBasis, weight::LaurentPolynomial)
    T = eltype(first(lpb))
    A = zeros(T, (length(lpb), length(lpb)))
    for I ∈ lpb.cross_index
        @inbounds A[I[1],I[2]] = scalar_product(weight * lpb[I[1]], lpb[I[2]], first(lpb.mesh), last(lpb.mesh))
        @inbounds A[I[2],I[1]] = A[I[1],I[2]] 
    end
    A
end

"""
    weight_mass_matrix(lpb::LaurentPolynomialBasis, , n::Int)

Special case of weighted mass matrix where the weight is the monomial X^n.
"""
function weight_mass_matrix(lpb::LaurentPolynomialBasis, n::Int)
    weight_mass_matrix(lpb::LaurentPolynomialBasis, Monomial(n))
end

"""
    deriv!(lpb::LaurentPolynomialBasis)

Deriv each elements of the basis.
"""
function deriv!(lpb::LaurentPolynomialBasis)
    for i in eachindex(lpb)
        deriv!(lpb[i])
    end
    lpb
end

"""
    deriv(lpb::LaurentPolynomialBasis)

Create an other basis made of the derivative of each elements of the basis.
"""
function deriv(lpb::LaurentPolynomialBasis)
    TL = eltype(lpb.elements)
    deriv_laurent = TL[]
    for i in eachindex(lpb)
        push!(deriv_laurent, deriv(lpb[i]))
    end
    LaurentPolynomialBasis(deriv_laurent)
end

"""
    build_on_basis(basis::LaurentPolynomialBasis{TL}, coeff) where TL

Compute the TL type elements (laurent Polynomial or piecewise laurent polynomial) on the basis
with respect to the coefficents.
"""
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