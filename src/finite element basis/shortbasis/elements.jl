########################################################################################
#                                  Default Elements
########################################################################################
struct DefaultElements{T} <: AbstractShortElements{T}
    polynomials::Vector{LaurentPolynomial{T}}
    size::Int
    binf::T
    bsup::T
    function DefaultElements(N, T, polynomials, size, binf, bsup)
        new{N,T}(polynomials, size, binf, bsup,)
    end
end

@inline Base.firstindex(::DefaultElements) = 1
@inline Base.lastindex(delem::DefaultElements) = delem.size
@inline Base.eachindex(delem::DefaultElements) = eachindex(delem.polynomials)
@inline Base.getindex(delem::DefaultElements, n::Int) =  delem.polynomials[n] 
@inline Base.first(delem::DefaultElements) = delem.polynomials[firstindex(delem)]
@inline Base.last(delem::DefaultElements) = delem.polynomials[lastindex(delem)]
@inline getpolynomial(delem::DefaultElements) = delem.polynomials

########################################################################################
#                                   P1 Elements
########################################################################################

struct P1Elements{T} <: AbstractShortElements{T}
    hfup::LaurentPolynomial{T}
    hfdown::LaurentPolynomial{T}
    dhfup::T
    dhfdown::T
    size::Int
    left::Bool
    right::Bool
    binf::T
    bsup::T
    function P1Elements(T::Type = Float64; left::Bool = false, right::Bool = false, binf::Real = -1, bsup::Real = 1)
        hfup = LaurentPolynomial([one(T),one(T)], 0, false, T(0))
        hfdown = LaurentPolynomial([one(T),-one(T)], 0, false, T(0))
        dhfup   = -one(T)
        dhfdown = one(T)
        new{T}(hfup, hfdown, dhfup, dhfdown, 2, left, right, binf, bsup)
    end
end

function ShortP1Basis(mesh::Mesh, T::Type = Float64;  kwargs...)
    p1elements = P1Elements(T; kwargs...)
    size = length(mesh) - 2 + p1elements.left + p1elements.right
    binf = p1elements.binf
    bsup = p1elements.bsup
    infos = InfoElement{T}[]
    if p1elements.left
        index = [2]
        m₁ = T(mesh[1])
        m₂ = T(mesh[2])
        segments = [1]
        shifts = [shift(T, m₁, m₂, binf, bsup)]
        invshifts = [shift(T, binf, bsup, m₁, m₂)]
        info = InfoElement(index, segments, shifts, invshifts)
        push!(infos, info)
    end
    for i ∈ eachindex(mesh)[2:end-1]
        index = [1,2]
        mᵢ₋₁ = T(mesh[i-1])
        mᵢ   = T(mesh[i])
        mᵢ₊₁ = T(mesh[i+1])
        segments = [i-1, i]
        shifts = [shift(T, mᵢ₋₁, mᵢ, binf, bsup), shift(T, mᵢ, mᵢ₊₁, binf, bsup)]
        invshifts = [shift(T, binf, bsup, mᵢ₋₁, mᵢ), shift(T, binf, bsup, mᵢ, mᵢ₊₁)]
        info = InfoElement(index, segments, shifts, invshifts)
        push!(infos, info)
    end
    if p1elements.right
        index = [1]
        mₑₙ₋₁ = T(mesh[end-1])
        mₑₙ = T(mesh[end])
        segments = [length(mesh)-1]
        shifts = [shift(T, mₑₙ₋₁, mₑₙ, binf, bsup)]
        invshifts = [shift(T, binf, bsup, mₑₙ₋₁, mₑₙ)]
        info = InfoElement(index, segments, shifts, invshifts)
        push!(infos, info)
    end

    ShortPolynomialBasis(p1elements, mesh, size, infos, 0)
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
@inline getpolynomial(p1::P1Elements) = [p1.hfup,p1.hfdown]
@inline getderivpolynomial(p1::P1Elements) = [p1.dhfup,p1.dhfdown]
########################################################################################
#                                  Integrated Legendre Elements
########################################################################################

struct IntLegendreElements{T} <: AbstractShortElements{T}
    polynomials::Vector{LaurentPolynomial{T}}
    derivpolynomials::Vector{LaurentPolynomial{T}}
    size::Int
    ordermin::Int
    ordermax::Int
    binf::T
    bsup::T

    function IntLegendreElements(T::Type = Float64; ordermin::Int = 2, ordermax = 2, binf::Real = -T(1), bsup::Real = T(1))
        @assert ordermin ≥ 1
        polynomials = LaurentPolynomial{T}[]
        derivpolynomials = LaurentPolynomial{T}[]  
        for n ∈ ordermin:ordermax
            Pₙ = Legendre(n-1; T = T, a = T(binf), b = T(bsup))
            push!(derivpolynomials, Pₙ)
            Qₙ = intLegendre(n-1; T = T, a = T(binf), b = T(bsup))
            push!(polynomials, Qₙ)
        end
        new{T}(polynomials, derivpolynomials, ordermax - ordermin + 1, ordermin, ordermax, T(binf), T(bsup))
    end
end

@inline Base.firstindex(::IntLegendreElements) = 1
@inline Base.eachindex(ilb::IntLegendreElements) = eachindex(ilb.polynomials)
@inline Base.getindex(ilb::IntLegendreElements, n::Int) =  ilb.polynomials[n] 
@inline getpolynomial(ilb::IntLegendreElements) = ilb.polynomials
@inline getderivpolynomial(ilb::IntLegendreElements) = ilb.derivpolynomials

function ShortIntLegendreBasis(mesh::Mesh, T::Type = Float64; Rcut::Real = last(mesh), normalize::Bool = false, kwargs...)
    intlegelement = IntLegendreElements(T; kwargs...)
    Ncut = min(findindex(mesh, Rcut), lastindex(mesh))
    size = intlegelement.size * (Ncut - 1)
    infos = Vector{InfoElement{T}}(undef, size)
    for i ∈ eachindex(mesh)[1:Ncut-1]
        segments = [i]
        shifts = [shift(T, mesh[i], mesh[i+1], intlegelement.binf, intlegelement.bsup)]
        invshifts = [shift(T, intlegelement.binf, intlegelement.bsup, mesh[i], mesh[i+1])]
        for n ∈ 1:intlegelement.size
            index = [n]
            infos[(i-1) * intlegelement.size + n] = InfoElement(index, segments, shifts, invshifts)
        end
    end
    ShortPolynomialBasis(intlegelement, mesh, size, infos, 0)
end

########################################################################################
#                                Difference Legendre Elements
# https://repositorio.ufba.br/bitstream/ri/13632/1/Marc%C3%ADlio%20N%20Guimar%C3%A3es.pdf
########################################################################################

struct DiffLegendreElements{T} <: AbstractShortElements{T}
    polynomials::Vector{LaurentPolynomial{T}}
    size::Int
    ordermax::Int
    binf::T
    bsup::T

    function DiffLegendreElements(T::Type = Float64; ordermax = 2, binf::Real = -1, bsup::Real = 1)
        polynomials = LaurentPolynomial{T}[]  
        for n ∈ 1:ordermax
            Qₙ = Legendre(n+1; T = T, a = binf, b = bsup) - Legendre(n-1; T = T, a = binf, b = bsup)
            push!(polynomials, Qₙ)
        end
        new{T}(polynomials, ordermax, ordermax, binf, bsup)
    end
end

@inline Base.firstindex(::DiffLegendreElements) = 1
@inline Base.eachindex(dlb::DiffLegendreElements) = eachindex(dlb.polynomials)
@inline Base.getindex(dlb::DiffLegendreElements, n::Int) =  dlb.polynomials[n] 
@inline getpolynomial(dlb::DiffLegendreElements) = dlb.polynomials

function ShortDiffLegendreBasis(mesh::Mesh, T::Type = Float64; normalize::Bool = false, kwargs...)
    difflegelement = DiffLegendreElements(T; kwargs...)
    size = difflegelement.size * (length(mesh) - 1)
    infos = Vector{InfoElement{T}}(undef, size)
    for i ∈ eachindex(mesh)[1:end-1]
        segments = [i]
        shifts = [shift(T, mesh[i], mesh[i+1], difflegelement.binf, difflegelement.bsup)]
        invshifts = [shift(T, difflegelement.binf, difflegelement.bsup, mesh[i], mesh[i+1])]
        for n ∈ 1:difflegelement.size
            index = [n]
            infos[(i-1) * difflegelement.size + n] = InfoElement(index, segments, shifts, invshifts)
        end
    end
    ShortPolynomialBasis(difflegelement, mesh, size, infos, 0)
end

