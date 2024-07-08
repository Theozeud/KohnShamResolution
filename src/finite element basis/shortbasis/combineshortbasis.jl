########################################################################################
#                            Combination Reduced Polynomial Basis
########################################################################################
struct InfoBlock
    index::CartesianIndex{2}
    rangerow::UnitRange{Int64}
    rangecolumn::UnitRange{Int64}
    diagonal::Bool
    interaction_index::Vector{CartesianIndex{2}}
end

@inline getindex(infoblock::InfoBlock) = infoblock.index
@inline getindex(infoblock::InfoBlock, i::Int) = infoblock.index[i]
@inline getrangerow(infoblock::InfoBlock) = infoblock.rangerow
@inline getrangecolumn(infoblock::InfoBlock) = infoblock.rangecolumn
@inline isdiagonal(infoblock::InfoBlock) = infoblock.diagonal
@inline getinteraction(infoblock::InfoBlock) = infoblock.interaction_index

struct CombineShortPolynomialBasis <: Basis
    basisVector
    size::Int
    blocks::Vector{InfoBlock}
    cumul_index::Vector{Int}
    function CombineShortPolynomialBasis(basisVector, size, blocks, cumul_index)
        new(basisVector, size, blocks, cumul_index)
    end
    function CombineShortPolynomialBasis(basisVector...)
        size = sum([length(basis) for basis ∈ basisVector])
        blocks = Vector{InfoBlock}(undef, length(basisVector)*(length(basisVector)+1)÷2)
        size_i = 1
        ib = 1
        cumul_index = zeros(Int, length(basisVector))
        for i ∈ eachindex(basisVector)
            cumul_index[i] = size_i
            size_j = 1
            for j ∈ 1:i
                rangerow = size_i:size_i+length(basisVector[i])-1
                rangecolumn = size_j:size_j+length(basisVector[j])-1
                size_j += length(basisVector[j])
                interaction_index = CartesianIndex{2}[]
                for n in eachindex(basisVector[i].infos)
                    for m in eachindex(basisVector[j].infos)
                        if !isempty(intersect(getsegments(basisVector[i], n), getsegments(basisVector[j], m)))
                            push!(interaction_index, CartesianIndex(n, m))
                        end
                    end
                end
                block = InfoBlock(CartesianIndex(i,j), rangerow, rangecolumn, i == j, interaction_index)
                blocks[ib] = block
                ib+=1
            end
            size_i += length(basisVector[i])
        end
        new(basisVector, size, blocks, cumul_index)
    end
end

@inline bottom_type(cb::CombineShortPolynomialBasis)= bottom_type(first(cb))

@inline Base.length(cb::CombineShortPolynomialBasis) = cb.size
@inline Base.eachindex(cb::CombineShortPolynomialBasis) = 1:cb.size
@inline Base.first(cb::CombineShortPolynomialBasis) = cb.basisVector[1]
@inline getbasis(cb::CombineShortPolynomialBasis, i::Int) = cb.basisVector[i]
@inline getblocks(cb::CombineShortPolynomialBasis) = cb.blocks

@inline getindex(cb::CombineShortPolynomialBasis, i::Int) =  getindex(cb.infoblock[i])
@inline getrangerow(cb::CombineShortPolynomialBasis, i::Int) = getrangerow(cb.infoblock[i])
@inline getrangecolumn(cb::CombineShortPolynomialBasis, i::Int) = getrangecolumn(cb.infoblock[i])
@inline isdiagonal(cb::CombineShortPolynomialBasis, i::Int) = isdiagonal(cb.infoblock[i])
@inline getinteraction(cb::CombineShortPolynomialBasis, i::Int) = getinteraction(cb.infoblock[i])

function find_basis(cb::CombineShortPolynomialBasis, i::Int)
    @assert i ≤ length(cb)
    ib = 1
    while (ib < length(cb.cumul_index)) && (cb.cumul_index[ib+1] ≤ i) 
        ib += 1
    end
    (ib, i-cb.cumul_index[ib]+1)
end


function mass_matrix(cb::CombineShortPolynomialBasis)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    for b ∈ getblocks(cb)
        @views ABlock = A[getrangerow(b), getrangecolumn(b)]
        if isdiagonal(b)
            fill_mass_matrix!(getbasis(cb, getindex(b,1)), ABlock)
        else
            fill_mass_matrix!(getbasis(cb, getindex(b,1)), getbasis(cb, getindex(b,2)), b.interaction_index, ABlock)
            @views ABlockT = A[getrangecolumn(b), getrangerow(b)]
            @. ABlockT = ABlock'
        end
    end
    A
end

function weight_mass_matrix(cb::CombineShortPolynomialBasis, weight::LaurentPolynomial)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    for b ∈ getblocks(cb)
        @views ABlock = A[getrangerow(b), getrangecolumn(b)]
        if isdiagonal(b)
            fill_weight_mass_matrix!(getbasis(cb, getindex(b,1)), weight, ABlock)
        else
            fill_weight_mass_matrix!(getbasis(cb, getindex(b,1)), getbasis(cb, getindex(b,2)), weight, b.interaction_index, ABlock)
            @views ABlockT = A[getrangecolumn(b), getrangerow(b)]
            @. ABlockT = ABlock'
        end
    end
    A
end

function weight_mass_matrix(cb::CombineShortPolynomialBasis, n::Int)
    weight_mass_matrix(cb, Monomial(n))
end

function fill_mass_matrix!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, interaction_index::Vector{CartesianIndex{2}}, A)
    for I ∈ interaction_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb1, I[1]), getsegments(spb2, I[2]))
            P = getpolynomial(spb1, I[1], i)
            Q = getpolynomial(spb2, I[2], j)
            invϕ = getinvshift(spb1, I[1], i)
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2]] += dinvϕ * scalar_product(P, Q, spb1.elements.binf, spb1.elements.bsup)
        end
        if isnormalized(spb1)
            @inbounds A[I[1], I[2]] *= getnormalization(spb1, I[1]) 
        end
        if isnormalized(spb2)
            @inbounds A[I[1], I[2]] *= getnormalization(spb2, I[2])
        end
    end
end

function fill_weight_mass_matrix!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, weight::LaurentPolynomial, interaction_index::Vector{CartesianIndex{2}}, A)
    for I ∈ interaction_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb1, I[1]), getsegments(spb2, I[2]))
            P = getpolynomial(spb1, I[1], i)
            Q = getpolynomial(spb2, I[2], j)
            invϕ = getinvshift(spb1, I[1], i)
            dinvϕ = invϕ[1]
            weight_shift = weight ∘ invϕ
            @inbounds A[I[1], I[2]] += dinvϕ * weight_scalar_product(P, Q, weight_shift, spb1.elements.binf, spb1.elements.bsup)
        end
        if isnormalized(spb1)
            @inbounds A[I[1], I[2]] *= getnormalization(spb1, I[1]) 
        end
        if isnormalized(spb2)
            @inbounds A[I[1], I[2]] *= getnormalization(spb2, I[2])
        end
    end
end

@memoize function build_basis(cb::CombineShortPolynomialBasis, i::Int)
    (ib, iib) = find_basis(cb, i)
    spb = getbasis(cb, ib)
    build_basis(spb, iib)
end

function build_on_basis(cb::CombineShortPolynomialBasis, coeffs)
    @assert eachindex(coeffs) == eachindex(cb) 
    poly = coeffs[firstindex(coeffs)] * build_basis(cb, firstindex(coeffs))
    for i ∈ eachindex(cb)[2:end]
        poly += coeffs[i] * build_basis(cb, i)
    end
    poly
end

function deriv(cb::CombineShortPolynomialBasis)
    derivBasisVector = [deriv(first(cb))]
    for i ∈ eachindex(cb.basisVector)[2:end]
        push!(derivBasisVector, deriv(getbasis(cb, i)))
    end
    CombineShortPolynomialBasis(derivBasisVector, cb.size, cb.blocks, cb.cumul_index)
end