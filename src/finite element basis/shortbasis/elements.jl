########################################################################################
#                                  Abstract Short Elements
########################################################################################
abstract type AbstractShortElements{N, T} end

@inline eltype(::AbstractShortElements{N, T}) where T = T
@inline isnormalized(::AbstractShortElements{N, T}) where {N,T} = N
@inline length(elem::AbstractShortElements) = elem.size
@inline getpolynomial(elem::AbstractShortElements, n::Int) = elem[n]


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
    if normalize
        normalization = T[]
        if left
            push!(normalization, √( 3 / (2 * (mesh[2] - mesh[1]) )  ) )
        end
        for i ∈ eachindex(mesh)[2:end-1]
            push!(normalization, √( 3 / (mesh[i+1] - mesh[i-1])  )    )
        end
        if right
            push!(normalization, √( 3 / (2 * (mesh[end] - mesh[end-1])))   )
        end
    else
        normalization = ones(size)
    end

    shifts = Vector{Shifts}[]
    binf = p1elements.binf
    bsup = p1elements.bsup
    if generate_shift
        if left
            push!(shifts, [generateshift(T, mesh[1], mesh[2], binf, bsup)])
        end
        for i ∈ eachindex(mesh)[2:end-1]
            shifts[i] = [generateshift(T, mesh[i-1], mesh[i], binf, bsup), generateshift(T, mesh[i], mesh[i+1], binf, bsup)]
        end
        if right
            push!(shifts, [generateshift(T, mesh[end-1], mesh[end], binf, bsup)])
        end
    else
        if left
            push!(shifts, [initiateshift(T, mesh[1], mesh[2], binf, bsup)])
        end
        for i ∈ eachindex(mesh)[2:end-1]
            shifts[i] = [initiateshift(T, mesh[i-1], mesh[i]), generateshift(T, mesh[i], mesh[i+1])]
        end
        if right
            push!(shifts, [initiateshift(T, mesh[end-1], mesh[end])])
        end
    end

    index = Int[]
    if left
        push!(index, [2])
    end
    for _ ∈ eachindex(mesh)[2:end-1]
        push!(index, [1,2])
    end
    if right
        push!(index, [1])
    end

    ShortPolynomialBasis(p1elements, mesh, size, normalization)
end

@inline Base.eachindex(::P1Elements) = 1:2
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

@inline Base.eachindex(ilb::IntLegendreElements) = eachindex(ilb.polynomials)
@inline Base.getindex(ilb::IntLegendreElements, n::Int) =  ilb.polynomials[n] 