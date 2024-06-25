mutable struct Shift{T}
    mᵢ::T
    mᵢ₊₁::T
    ϕ::Union{LaurentPolynomial{T}, Nothing}
end

function generateshift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    Shift(T(mᵢ), T(mᵢ₊₁), shift(T, a, b, mᵢ, mᵢ₊₁))
end
function initiateshift(T::Type, mᵢ::Real, mᵢ₊₁::Real)
    Shift(T(mᵢ), T(mᵢ₊₁), nothing)
end

struct ShortPolynomialBasis{TB} <: Basis
    elements::TB
    mesh::OneDMesh
    size::Int
    normalization
    shifts
    index::Vector{Int}
    infos_integrate
end

@inline eltype(::ShortPolynomialBasis{TB}) where TB = TB
@inline bottom_type(spb::ShortPolynomialBasis) = eltype(spb.elements)
@inline Base.length(spb::ShortPolynomialBasis) = spb.size

@inline function getshift(spb::ShortPolynomialBasis, i::Int, j::Int)
    if spb.shifts[i][j].ϕ isa Nothing
        generateshift(spb, i, j)
    end
    spb.shifts[i][j].ϕ
end

@inline function generateshift(spb, i, j)
    T = bottom_type(spb)
    spb.shifts[i][j].ϕ = shift(T, spb.basis.binf, spb.basis.bsup, spb.shifts[i].mᵢ, spb.shifts[i].mᵢ₊₁)
    nothing
end

## How to define this functions properly ?
#@inline Base.getindex(spb::ShortPolynomialBasis, n::Int) =  spb.elements[n] 
#@inline Base.iterate(spb::ShortPolynomialBasis, state = 1) = state > length(spb) ? nothing : (spb[state], state+1)
@inline Base.eachindex(spb::ShortPolynomialBasis) = 1:spb.size

function mass_matrix(spb::ShortPolynomialBasis)
    @unpack elements, mesh, size = spb
    T = bottom_type(spb)
    A = zeros(T, size, size)
    fill_mass_matrix!(spb, A)
    A
end

function fill_mass_matrix!(spb::ShortPolynomialBasis, A)
    @unpack elements, mesh, infos_integrate = spb
    fill_mass_matrix!(elements, mesh, A)
    nothing
end

function weight_mass_matrix(spb::ShortPolynomialBasis, weight::LaurentPolynomial)
    T = bottom_type(spb)
    A = zeros(T, spb.size, spb.size)
    fill_weight_mass_matrix!(spb, A, weight)
    A
end

function fill_weight_mass_matrix!(spb::ShortPolynomialBasis, A, weight::LaurentPolynomial)
    @unpack elements, normalization = spb
    for I ∈ filled_index
        for (ip,iq, ϕ)∈ getinfo_integrate(spb, I)
            P = getpolynomial(elements, ip)
            Q = getpolynomial(elements, iq)
            weight_shift = weight ∘ ϕ
            @inbounds A[I[1], I[2]] += weight_scalar_product(P, Q, weight_shift, basis.binf, basis.bsup)
        end
        @inbounds A[I[1], I[2]] *= normalization(I[1]) * normalization(I[2])
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

function build_on_basis(spb::ShortPolynomialBasis, coeffs)
    @assert eachindex(coeffs) == eachindex(spb) 
    poly = zero(first(spb.elements)) #???
    for i ∈ eachindex(spb)
        for j ∈ eachindex(index[i])
            ϕ = getshift(spb, i, j)
            poly += coeff[i] * spb.normalization[i] * getpolynomial(spb.elements, j) ∘ ϕ
        end
    end
    poly
end
