########################################################################################
#                            Combination Reduced Polynomial Basis
########################################################################################
struct CombineShortPolynomialBasis{T}
    basisVector::Vector{ShortPolynomialBasis}
    blocks::Vector{Tuple{CartesianIndex{2}, UnitRange{Int64}, UnitRange{Int64}, Bool}}
    size::Int
    interaction::Bool
end

function CombineBasis(basisVector...; interaction = true)
    blocks = Vector{Tuple{CartesianIndex{2}, UnitRange{Int64}, UnitRange{Int64}, Bool}}[]
    size_i = 1
    for i ∈ eachindex(basisVector)
        if interaction
            size_j = 1
            for j ∈ 1:i
                irange = size_i:size_i+size(basisVector[i])
                jrange = size_j:size_j+size(basisVector[j])
                push!(blocks, (CartesianIndex(i,j), irange, jrange, j==i))
                size_j += size(basisVector[j])
            end
        else
            irange = size_i:size_i+size(basisVector[i])
            push!(blocks, (CartesianIndex(i,j), irange, irange, true))
        end
        size_i += size(basisVector[i])
    end
    CombineShortPolynomialBasis(basisVector, blocks, size_i, interaction)
end

@inline Base.length(cb::CombineShortPolynomialBasis) = cb.size
@inline bottom_type(cb::CombineShortPolynomialBasis) = eltype(first(cb.polynomials))
@inline is_diagonal_block(cb::CombineShortPolynomialBasis, ib) = cb.blocks[ib][4] == true
@inline getbasis(cb::CombineShortPolynomialBasis, i::Int) = cb.basis[i]

function mass_matrix(cb::CombineShortPolynomialBasis)
    T = eltype(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    for (ib,b) ∈ enumerate(cb.blocks)
        @views ABlock = A[b[2], b[3]]
        if is_diagonal_block(rb, ib)
            fill_mass_matrix!(getbasis(cb, b[1][1]), ABlock)
        else
            fill_mass_matrix!(getbasis(cb, b[1][1]), getbasis(cb, b[1][2]), ABlock)
            @views ABlockT = A[b[3], b[2]]
            @. ABlockT = ABlock
        end
    end
    A
end

function weight_mass_matrix(cb::CombineShortPolynomialBasis, weight::LaurentPolynomial)
    T = eltype(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    for (ib,b) ∈ enumerate(cb.blocks)
        @views ABlock = A[b[2], b[3]]
        if is_symmetrical_block(rb, ib)
            fill_weight_mass_matrix!(getbasis(cb, b[1][1]), ABlock, weight)
        else
            fill_weight_mass_matrix!(getbasis(cb, b[1][1]), getbasis(cb, b[1][2]), ABlock, weight)
            @views ABlockT = A[b[2], b[1]]
            @. ABlockT = ABlock
        end
    end
    A
end