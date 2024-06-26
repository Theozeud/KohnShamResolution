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
@inline Base.iterate(ielem::InfoElement, state = 1) = state > length(ielem) ? nothing : ((getindex(ielem, i), getsegments(ielem, i), getshift(ielem, i), getinvshift(ielem, i)), state+1)

@inline arecompatible(ielem1::InfoElement, ielem2::InfoElement) = !isempty(intersect(ielem1.segments, ielem2.segments))

struct ShortPolynomialBasis{TB} <: Basis
    elements::TB
    mesh::OneDMesh
    size::Int
    infos::Vector{InfoElement}
    coupling_index::CartesianIndex{2}

    function ShortPolynomialBasis(elements, mesh::OneDMesh, size::Int, infos::Vector{InfoElement})
        
        # Basics checks for consistency
        @assert size == length(infos)
        for info ∈ infos
            @assert getsegments(info) ⊆ eachindex(mesh)[1:end-1]
        end 

        # Creating the coupling indices
        coupling_index = CartesianIndex{2}[]
        for i in eachindex(infos)
            for j in 1:i
                if !isempty(intersect(getsegment(infos[i]), getsegment(infos[j])))
                    push!(coupling_index, CartesianIndex(i, j))
                end
            end
        end
        new{typeof(elements)}(elements, mesh, size, infos, coupling_index)
    end
end

@inline eltype(::ShortPolynomialBasis{TB}) where TB = TB
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
    for I ∈ spb.coupling_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb, I[1]), getsegments(spb, I[2]))
            P = getpolynomial(spb, I[1], i)
            Q = getpolynomial(spb, I[2], j)
            ϕ = getshift(spb, I[1], i)
            weight_shift = weight ∘ ϕ
            @inbounds A[I[1], I[2]] += weight_scalar_product(P, Q, weight_shift, basis.binf, basis.bsup)
        end
        @inbounds A[I[1], I[2]] *= getnormalization(spb, I[1]) * getnormalization(spb, I[2])
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

function build_on_basis(spb::ShortPolynomialBasis, coeffs)
    @assert eachindex(coeffs) == eachindex(spb) 
    poly = zero(first(spb.elements))
    for i ∈ eachindex(spb)
        for (j,_,_,invϕ) ∈ eachindex(spb.infos[i])
            poly += coeff[i] * getnormalization(spb,i) * getpolynomial(spb.elements, j) ∘ invϕ
        end
    end
    poly
end
