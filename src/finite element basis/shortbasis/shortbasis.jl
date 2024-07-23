########################################################################################
#                                  Info Element
########################################################################################
struct InfoElement{T}
    index::Vector{Int}
    segments::Vector{Int}
    ϕ::Vector{LaurentPolynomial{T}}
    invϕ::Vector{LaurentPolynomial{T}}
end

@inline _getindex(ielem::InfoElement) = ielem.index
@inline _getindex(ielem::InfoElement, i::Int) = ielem.index[i]
@inline getsegments(ielem::InfoElement) = ielem.segments
@inline getsegments(ielem::InfoElement, i::Int) = ielem.segments[i]
@inline getshift(ielem::InfoElement, i::Int) = ielem.ϕ[i]
@inline getinvshift(ielem::InfoElement, i::Int) = ielem.invϕ[i]
@inline getnormalization(ielem::InfoElement, i::Int) = ielem.invϕ[i][1]

@inline Base.length(ielem::InfoElement) = length(ielem.index)
@inline Base.eachindex(ielem::InfoElement) = eachindex(ielem.index)
@inline Base.iterate(ielem::InfoElement, state = 1) = state > length(ielem) ? nothing : ((_getindex(ielem, state), getsegments(ielem, state), getshift(ielem, state), getinvshift(ielem, state)), state+1)

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
    deriv_order::Int

    function ShortPolynomialBasis(elements, mesh::OneDMesh, size::Int, infos, coupling_index::Vector{CartesianIndex{2}}, deriv_order::Int) 
        new{typeof(elements)}(elements, mesh, size, infos, coupling_index, deriv_order)
    end

    function ShortPolynomialBasis(elements, mesh::OneDMesh, size::Int, infos, deriv_order) 
        
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
        new{typeof(elements)}(elements, mesh, size, infos, coupling_index, deriv_order)
    end
end

@inline Base.eltype(::ShortPolynomialBasis{TB}) where TB = TB
@inline bottom_type(spb::ShortPolynomialBasis) = eltype(spb.elements)
@inline Base.length(spb::ShortPolynomialBasis) = spb.size
@inline Base.eachindex(spb::ShortPolynomialBasis) = 1:spb.size

@inline isnormalized(spb::ShortPolynomialBasis) = isnormalized(spb.elements)

@inline _getindex(spb::ShortPolynomialBasis, i::Int) = _getindex(spb.infos[i])
@inline _getindex(spb::ShortPolynomialBasis, i::Int, j::Int) = _getindex(spb.infos[i], j)
@inline getsegments(spb::ShortPolynomialBasis, i::Int) = getsegments(spb.infos[i])
@inline getsegments(spb::ShortPolynomialBasis, i::Int, j::Int) = getsegments(spb.infos[i], j)
@inline getshift(spb::ShortPolynomialBasis, i::Int, j::Int) = getshift(spb.infos[i], j)
@inline getinvshift(spb::ShortPolynomialBasis, i::Int, j::Int) = getinvshift(spb.infos[i], j)

@inline getbasis(spb::ShortPolynomialBasis, ::Int) = spb
function find_basis(spb::ShortPolynomialBasis, i::Int)
    @assert i ≤ length(spb)
    (1,i)
end

@inline function getnormalization(spb::ShortPolynomialBasis, i::Int)
    norma = bottom_type(spb)(0)
    for j ∈ eachindex(spb.infos[i])
        norma += getnormalization(spb.infos[i], j)^(1+spb.deriv_order*2) * getnormalization(spb.elements, _getindex(spb, i, j))
    end
    1/sqrt(norma)
end

@inline getpolynomial(spb::ShortPolynomialBasis, i::Int, j::Int)= getpolynomial(spb.elements, _getindex(spb, i, j))


########################################################################################
#                                     Mass Matrix && Co
########################################################################################

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
            if isnormalized(spb)
                @inbounds A[I[1], I[2]] *= getnormalization(spb, I[1]) * getnormalization(spb, I[2])
            end
            @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
        end
    end
    nothing
end

function weight_mass_matrix(spb::ShortPolynomialBasis, weight)
    T = bottom_type(spb)
    A = zeros(T, spb.size, spb.size)
    fill_weight_mass_matrix!(spb, weight, A)
    A
end

function weight_mass_matrix(spb::ShortPolynomialBasis, n::Int)
    weight_mass_matrix(spb, Monomial(n))
end

function fill_weight_mass_matrix!(spb::ShortPolynomialBasis, weight, A)
    for I ∈ spb.coupling_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb, I[1]), getsegments(spb, I[2]))
            P = getpolynomial(spb, I[1], i)
            Q = getpolynomial(spb, I[2], j)
            invϕ = getinvshift(spb, I[1], i)
            weight_shift = weight ∘ invϕ
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2]] += dinvϕ * weight_scalar_product(P, Q, weight_shift, spb.elements.binf, spb.elements.bsup)
        end
        if isnormalized(spb)
            @inbounds A[I[1], I[2]] *= getnormalization(spb, I[1]) * getnormalization(spb, I[2])
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

function weight_mass_vector(spb::ShortPolynomialBasis, weight)
    T = bottom_type(spb)
    A = zeros(T, spb.size)
    fill_weight_mass_vector!(spb, weight, A)
    A
end

function fill_weight_mass_vector!(spb::ShortPolynomialBasis, weight, A)
    for i ∈ eachindex(spb)
        for j ∈ eachindex(spb.infos[i])
            P = getpolynomial(spb, i, j)
            invϕ = getinvshift(spb, i, j)
            weight_shift = weight ∘ invϕ
            dinvϕ = invϕ[1]
            @inbounds A[i] += dinvϕ * weight_scalar_product(P, weight_shift, spb.elements.binf, spb.elements.bsup)
        end
        if isnormalized(spb)
            @inbounds A[i] *= getnormalization(spb, i)
        end
    end
end

@memoize function build_basis(spb::ShortPolynomialBasis, i::Int)
    T = bottom_type(spb)
    polys = LaurentPolynomial{T}[]
    for (j,_,ϕ,_) ∈ spb.infos[i]
        if isnormalized(spb)
            push!(polys, getnormalization(spb,i) * getpolynomial(spb.elements, j) ∘ ϕ)
        else
            push!(polys, getpolynomial(spb.elements, j) ∘ ϕ)
        end
    end
    PiecewiseLaurentPolynomial(spb.mesh, polys, getsegments(spb, i), T(0))
end

function build_on_basis(spb::ShortPolynomialBasis, coeffs)
    @assert eachindex(coeffs) == eachindex(spb) 
    poly = coeffs[firstindex(coeffs)] * build_basis(spb, firstindex(coeffs))
    for i ∈ eachindex(spb)[2:end]
        poly += coeffs[i] * build_basis(spb, i)
    end
    poly
end

function deriv(spb::ShortPolynomialBasis)
    deriv_polynomials = deriv.(getpolynomial(spb.elements))
    deriv_elements = DefaultElements(isnormalized(spb.elements), eltype(spb.elements), deriv_polynomials, spb.elements.size, spb.elements.binf, spb.elements.bsup, spb.elements.normalization)
    ShortPolynomialBasis(deriv_elements, spb.mesh, spb.size, spb.infos, spb.coupling_index, spb.deriv_order+1)
end