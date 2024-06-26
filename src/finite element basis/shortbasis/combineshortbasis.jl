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

struct CombineShortPolynomialBasis
    basisVector::Vector{ShortPolynomialBasis}
    size::Int
    blocks::Vector{InfoBlock}
    function CombineBasis(basisVector...)
        size = length([length(basis) for basis ∈ basisVector])
        blocks = Vector{InfloBlock}(undef, length(basisVector))
        size_i = 1
        for i ∈ eachindex(basisVector)
            size_j = 1
            for j ∈ 1:i
                rangerow = size_i:size_i+size(basisVector[i])
                rangecolumn = size_j:size_j+size(basisVector[j])
                size_j += size(basisVector[j])
                interaction_index = CartesianIndex{2}[]
                for n in eachindex(basisVector[i].infos)
                    for m in eachindex(basisVector[j].infos)
                        if !isempty(intersect(getsegment(basisVector[i].infos, n), getsegment(basisVector[j].infos[m])))
                            push!(coupling_index, CartesianIndex(n, m))
                        end
                    end
                end
                block = InfoBlock(CartesianIndex(i,j), rangerow, rangecolumn, i == j, interaction_index)
                blocks[i] = block
            end
        end
        new(basisVector, size, blocks)
    end
end

@inline getbasis(cb::CombineShortPolynomialBasis, i::Int) = cb.basis[i]
@inline Base.length(cb::CombineShortPolynomialBasis) = cb.size
@inline getblocks(cb::CombineShortPolynomialBasis) = cb.block

@inline getindex(cb::CombineShortPolynomialBasis, i::Int) =  getindex(cb.infoblock[i])
@inline getrangerow(cb::CombineShortPolynomialBasis, i::Int) = getrangerow(cb.infoblock[i])
@inline getrangecolumn(cb::CombineShortPolynomialBasis, i::Int) = getrangecolumn(cb.infoblock[i])
@inline isdiagonal(cb::CombineShortPolynomialBasis, i::Int) = isdiagonal(cb.infoblock[i])
@inline getinteraction(cb::CombineShortPolynomialBasis, i::Int) = getinteraction(cb.infoblock[i])

function mass_matrix(cb::CombineShortPolynomialBasis)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    for b ∈ getblocks(cb)
        @views ABlock = A[getrangerow(b), getrangecolumn(b)]
        if isdiagonal(b)
            fill_mass_matrix!(getbasis(cb, getindex(b,1)), ABlock)
        else
            fill_mass_matrix!(getbasis(cb, getindex(b,1)), getbasis(cb, getindex(b,2)), cb.interaction_index, A)
            @views ABlockT = A[getrangecolumn(b), getrangerow(b)]
            @. ABlockT = ABlock
        end
    end
    A
end

function weight_mass_matrix(cb::CombineShortPolynomialBasis, weight::LaurentPolynomial)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    for b ∈ enumerate(cb.blocks)
        @views ABlock = A[getrangerow(b), getrangecolumn(b)]
        if isdiagonal(b)
            fill_weight_mass_matrix!(getbasis(cb, getindex(b,1)), weight, ABlock)
        else
            fill_weight_mass_matrix!(getbasis(cb, getindex(b,1)), getbasis(cb, getindex(b,2)), weight, cb.interaction_index, A)
            @views ABlockT = A[getrangecolumn(b), getrangerow(b)]
            @. ABlockT = ABlock
        end
    end
    A
end

function fill_mass_matrix!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, interaction_index::Vector{CartesianIndex{2}}, A)
    for I ∈ interaction_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb1, I[1]), getsegments(spb2, I[2]))
            P = getpolynomial(spb1, I[1], i)
            Q = getpolynomial(spb2, I[2], j)
            @inbounds A[I[1], I[2]] += scalar_product(P, Q, basis.binf, basis.bsup)
        end
        @inbounds A[I[1], I[2]] *= getnormalization(spb1, I[1]) * getnormalization(spb2, I[2])
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
end

function fill_weight_mass_matrix!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, weight::LaurentPolynomial, interaction_index::Vector{CartesianIndex{2}}, A)
    for I ∈ interaction_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb1, I[1]), getsegments(spb2, I[2]))
            P = getpolynomial(spb1, I[1], i)
            Q = getpolynomial(spb2, I[2], j)
            ϕ = getshift(spb1, I[1], i)
            weight_shift = weight ∘ ϕ
            @inbounds A[I[1], I[2]] += weight_scalar_product(P, Q, weight_shift, basis.binf, basis.bsup)
        end
        @inbounds A[I[1], I[2]] *= getnormalization(spb1, I[1]) * getnormalization(spb2, I[2])
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
end