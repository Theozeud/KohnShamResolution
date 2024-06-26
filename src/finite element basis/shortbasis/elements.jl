########################################################################################
#                                  Abstract Short Elements
########################################################################################
abstract type AbstractShortElements{N, T} end

@inline Base.eltype(::AbstractShortElements{N, T}) where T = T
@inline isnormalized(::AbstractShortElements{N, T}) where {N,T} = N
@inline Base.length(elem::AbstractShortElements) = elem.size
@inline getpolynomial(elem::AbstractShortElements, n::Int) = elem[n]
@inline Base.first(elem::AbstractShortElements) = elem[firstindex(elem)]

########################################################################################
#                                   P1 Elements
########################################################################################

struct P1Elements{N, T} <: AbstractShortElements{N, T}
    hfup::LaurentPolynomial{T}
    hfdown::LaurentPolynomial{T}
    size::Int
    left::Bool
    right::Bool
    binf::T
    bsup::T
    function P1Elements(T::Type = Float64; left::Bool = false, right::Bool = false, normalize::Bool = false, binf::Real = -1, bsup::Real = 1)
        hfup = LaurentPolynomial([T(1)], 1, false, T(0))
        hfdown = LaurentPolynomial([T(-1)], 1, false, T(0))
        new{normalize, T}(hfup, hfdown, 2, left, right, binf, bsup)
    end
end

function ShortP1Basis(mesh::OneDMesh, T::Type = Float64; normalize::Bool = false, generate_shift = true, kwargs...)
    p1elements = P1Elements(T; normalize = normalize, kwargs)
    size = length(mesh) - 2 + left + right
    binf = p1elements.binf
    bsup = p1elements.bsup
    infos = InfoElement{T}[]
    if left
        index = [2]
        m₁ = T(mesh[2])
        m₂ = T(mesh[1])
        normalization = normalize ? √( T(3) / (T(2) * (m₂ - m₁))) : T(1)
        segments = [1]
        shifts = [shift(T, m₁, m₂, binf, bsup)]
        info = InfoElement(index, segments, shifts, normalization)
        push!(infos, info)
    end
    for i ∈ eachindex(mesh)[2:end-1]
        index = [1,2]
        mᵢ₋₁ = T(mesh[i-1])
        mᵢ   = T(mesh[i])
        mᵢ₊₁ = T(mesh[i+1])
        normalization = normalize ? √( T(3) / (mᵢ₊₁ - mᵢ₋₁) ) : T(1)
        segments = [i-1, i]
        shifts = [shift(T, mᵢ₋₁, mᵢ, binf, bsup), shift(T, mᵢ, mᵢ₊₁, binf, bsup)]
        info = InfoElement(index, segments, shifts, normalization)
        push!(infos, info)
    end
    if right
        index = [1]
        mₑₙ₋₁ = T(mesh[end-1])
        mₑₙ = T(mesh[end])
        normalization = normalize ? √( 3 / (2 * (mₑₙ - mₑₙ₋₁))) : T(1)
        segments = [length(mesh)-1]
        shifts = [shift(T, mₑₙ₋₁, mₑₙ, binf, bsup)]
        info = InfoElement(index, segments, shifts, normalization)
        push!(infos, info)
    end

    ShortPolynomialBasis(p1elements, mesh, size, infos)
end

@inline Base.eachindex(::P1Elements) = 1:2
@inline Base.firstindex(::P1Elements) = 1
@inline function Base.getindex(p1::P1Elements, n::Int)
    if n == 1
        return p1.hfup
    elseif n == 2
        return p1.hfdown
    else
        throw(BoundsError(p1, n))
    end
end

########################################################################################
#                                  Integrated Legendre Elements
########################################################################################

struct IntLegendreElements{N, T} <: AbstractShortElements{N, T}
    polynomials::Vector{LaurentPolynomial{T}}
    size::Int
    ordermin::Int
    ordermax::Int
    binf::T
    bsup::T
    interaction

    function IntLegendreElements(T::Type = Float64; ordermin::Int = 2, ordermax = 2, normalize::Bool = false, binf::Real = -1, bsup::Real = 1)
        @assert ordermin ≥ 1
        polynomials = LaurentPolynomial{T}[]    
        for n ∈ ordermin:ordermax
            Qₙ = intLegendre(n-1; T = T, normalize = true, a = binf, b = bsup)
            push!(polynomials, Qₙ)
        end
        new{normalize, T}(polynomials, ordermax - ordermin + 1, ordermin, ordermax, binf, bsup, Matrix(I, ordermax - ordermin + 1, ordermax - ordermin + 1))
    end
end

@inline Base.firstindex(::IntLegendreElements) = 1
@inline Base.eachindex(ilb::IntLegendreElements) = eachindex(ilb.polynomials)
@inline Base.getindex(ilb::IntLegendreElements, n::Int) =  ilb.polynomials[n] 


function ShortIntLegendreBasis(mesh::OneDMesh, T::Type = Float64; normalize::Bool = false, generate_shift = true, kwargs...)
    intlegelement = IntLegendreElements(T; normalize = normalize, kwargs...)
    size = intlegelement.size * (length(mesh) - 1)
    infos = Vector{InfoElement{T}}(undef, size)
    for i ∈ eachindex(mesh)[1:end-1]
        normalization = normalize ? sqrt((intlegelement.bsup - intlegelement.binf)/ (mesh[i+1] - mesh[i])) : T(1)
        segments = [i]
        shifts = [shift(T, intlegelement.binf, intlegelement.bsup, mesh[i], mesh[i])]
        for n ∈ 1:intlegelement.size
            index = [n]
            infos[(i-1) * intlegelement.size + n] = InfoElement(index, segments, shifts, normalization)
        end
    end
    ShortPolynomialBasis(intlegelement, mesh, size, infos)
end