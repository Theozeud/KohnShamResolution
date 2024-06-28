########################################################################################
#                                  Info Element
########################################################################################
struct InfoElement{T}
    index::Vector{Int}
    segments::Vector{Int}
    ϕ::Vector{LaurentPolynomial{T}}
    invϕ::Vector{LaurentPolynomial{T}}
    normalization::T
end

@inline getindex(ielem::InfoElement) = ielem.index
@inline getindex(ielem::InfoElement, i::Int) = ielem.index[i]
@inline getsegments(ielem::InfoElement) = ielem.segments
@inline getsegments(ielem::InfoElement, i::Int) = ielem.segments[i]
@inline getshift(ielem::InfoElement, i::Int) = ielem.ϕ[i]
@inline getinvshift(ielem::InfoElement, i::Int) = ielem.invϕ[i]
@inline getnormalization(ielem::InfoElement) = ielem.normalization

@inline Base.length(ielem::InfoElement) = length(ielem.index)
@inline Base.iterate(ielem::InfoElement, state = 1) = state > length(ielem) ? nothing : ((getindex(ielem, state), getsegments(ielem, state), getshift(ielem, state), getinvshift(ielem, state)), state+1)

@inline arecompatible(ielem1::InfoElement, ielem2::InfoElement) = !isempty(intersect(ielem1.segments, ielem2.segments))

########################################################################################
#                                     Short Basis
########################################################################################
struct ShortPolynomialBasis{TB <: AbstractShortElements} <: Basis
    elements::TB
    mesh::OneDMesh
    size::Int
    infos::Vector{InfoElement}
    coupling_index::Vector{CartesianIndex{2}}

    function ShortPolynomialBasis(elements, mesh::OneDMesh, size::Int, infos, coupling_index::Vector{CartesianIndex{2}}) 
        new{typeof(elements)}(elements, mesh, size, infos, coupling_index)
    end

    function ShortPolynomialBasis(elements, mesh::OneDMesh, size::Int, infos) 
        
        # Basics checks for consistency
        @assert size == length(infos)
        for info ∈ infos
            @assert getsegments(info) ⊆ eachindex(mesh)[1:end-1]
        end 

        # Creating the coupling indices
        coupling_index = CartesianIndex{2}[]
        for i in eachindex(infos)
            for j in 1:i
                if !isempty(intersect(getsegments(infos[i]), getsegments(infos[j])))
                    push!(coupling_index, CartesianIndex(i, j))
                end
            end
        end
        new{typeof(elements)}(elements, mesh, size, infos, coupling_index)
    end
end

@inline Base.eltype(::ShortPolynomialBasis{TB}) where TB = TB
@inline bottom_type(spb::ShortPolynomialBasis) = eltype(spb.elements)
@inline Base.length(spb::ShortPolynomialBasis) = spb.size

@inline getindex(spb::ShortPolynomialBasis, i::Int) = getindex(spb.infos[i])
@inline getsegments(spb::ShortPolynomialBasis, i::Int) = getsegments(spb.infos[i])
@inline getindex(spb::ShortPolynomialBasis, i::Int, j::Int) = getindex(spb.infos[i], j)
@inline getsegments(spb::ShortPolynomialBasis, i::Int, j::Int) = getsegments(spb.infos[i], j)
@inline getshift(spb::ShortPolynomialBasis, i::Int, j::Int) = getshift(spb.infos[i], j)
@inline getinvshift(spb::ShortPolynomialBasis, i::Int, j::Int) = getinvshift(spb.infos[i], j)
@inline getnormalization(spb::ShortPolynomialBasis, i::Int) = getnormalization(spb.infos[i])

@inline getpolynomial(spb::ShortPolynomialBasis, i::Int, j::Int)= getpolynomial(spb.elements, getindex(spb, i, j))

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
    @unpack elements, mesh = spb
    try 
        fill_mass_matrix!(elements, mesh, A)
    catch
        for I ∈ spb.coupling_index
            for (i,j) ∈ intersection_with_indices(getsegments(spb, I[1]), getsegments(spb, I[2]))
                P = getpolynomial(spb, I[1], i)
                Q = getpolynomial(spb, I[2], j)
                invϕ = getinvshift(spb, I[1], i)
                dinvϕ = invϕ[1]
                @inbounds A[I[1], I[2]] += dinvϕ * scalar_product(P, Q, spb.elements.binf, spb.elements.bsup)
            end
            @inbounds A[I[1], I[2]] *= getnormalization(spb, I[1]) * getnormalization(spb, I[2])
            @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
        end
    end
    nothing
end

function weight_mass_matrix(spb::ShortPolynomialBasis, weight::LaurentPolynomial)
    T = bottom_type(spb)
    A = zeros(T, spb.size, spb.size)
    fill_weight_mass_matrix!(spb, weight, A)
    A
end

function weight_mass_matrix(spb::ShortPolynomialBasis, n::Int)
    weight_mass_matrix(spb, Monomial(n))
end

function fill_weight_mass_matrix!(spb::ShortPolynomialBasis, weight::LaurentPolynomial, A)
    for I ∈ spb.coupling_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb, I[1]), getsegments(spb, I[2]))
            P = getpolynomial(spb, I[1], i)
            Q = getpolynomial(spb, I[2], j)
            invϕ = getinvshift(spb, I[1], i)
            weight_shift = weight ∘ invϕ
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2]] += dinvϕ * weight_scalar_product(P, Q, weight_shift, spb.elements.binf, spb.elements.bsup)
        end
        @inbounds A[I[1], I[2]] *= getnormalization(spb, I[1]) * getnormalization(spb, I[2])
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

@memoize function build_basis(spb::ShortPolynomialBasis, i::Int)
    T = bottom_type(spb)
    polys = LaurentPolynomial{T}[]
    for (j,_,ϕ,_) ∈ spb.infos[i]
        push!(polys, getnormalization(spb,i) * getpolynomial(spb.elements, j) ∘ ϕ)
    end
    PiecewiseLaurentPolynomial(spb.mesh, polys, getsegments(spb, i), T(0))
end

function build_on_basis(spb::ShortPolynomialBasis, coeffs)
    @assert eachindex(coeffs) == eachindex(spb) 
    poly = coeff[firstindex(coeffs)] * build_basis(spb, firstindex(coeffs))
    for i ∈ eachindex(spb)[2:end]
        poly += coeff[i] * build_basis(spb, i)
    end
    poly
end

function deriv(spb::ShortPolynomialBasis)
    deriv_elements = DefaultElements(isnormalized(spb.elements), eltype(spb.elements), deriv.(getpolynomial(spb.elements)), spb.elements.size, spb.elements.binf, spb.elements.bsup)
    ShortPolynomialBasis(deriv_elements, spb.mesh, spb.size, spb.infos, spb.coupling_index)
end