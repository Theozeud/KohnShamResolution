########################################################################################
#                            Combination Reduced Polynomial Basis
########################################################################################
struct InfoBlock{N}
    index::CartesianIndex{N}
    axes::NTuple{N,UnitRange{Int64}}
    diagonal::Bool
    interaction_index::Vector{CartesianIndex{N}}
end

@inline _getindex(infoblock::InfoBlock) = infoblock.index
@inline _getindex(infoblock::InfoBlock, i::Int) = infoblock.index[i]
@inline getrangerow(infoblock::InfoBlock) = infoblock.axes[1]
@inline getrangecolumn(infoblock::InfoBlock) = infoblock.axes[2]
@inline getaxes(infoblock::InfoBlock, n::Int) = infoblock.axes[n]
@inline isdiagonal(infoblock::InfoBlock) = infoblock.diagonal
@inline getinteraction(infoblock::InfoBlock) = infoblock.interaction_index

struct CombineShortPolynomialBasis <: Basis
    basisVector
    size::Int
    blocks::Vector{InfoBlock}
    blocks3::Vector{InfoBlock}
    cumul_index::Vector{Int}
    function CombineShortPolynomialBasis(basisVector, size, blocks, blocks3, cumul_index)
        new(basisVector, size, blocks, blocks3, cumul_index)
    end
    function CombineShortPolynomialBasis(basisVector...)
        size = sum([length(basis) for basis ∈ basisVector])
        # Infos Block for the matrices
        blocks = Vector{InfoBlock}(undef, length(basisVector)*(length(basisVector)+1)÷2)
        blocks3 = Vector{InfoBlock}(undef, length(basisVector)*(length(basisVector)+1)*(length(basisVector)+2)÷6)
        size_i = 1
        ib = 1
        ib3 = 1
        cumul_index = zeros(Int, length(basisVector))
        for i ∈ eachindex(basisVector)
            cumul_index[i] = size_i
            axes1 = size_i:size_i+length(basisVector[i])-1
            size_j = 1
            for j ∈ 1:i
                axes2 = size_j:size_j+length(basisVector[j])-1
                interaction_index = CartesianIndex{2}[]
                size_k = 1
                for k ∈ 1:j
                    axes3 = size_k:size_k+length(basisVector[k])-1
                    interaction_index3 = CartesianIndex{3}[] 
                    for n in eachindex(basisVector[i].infos)
                        for m in eachindex(basisVector[j].infos)
                            intersect_nm = intersect(getsegments(basisVector[i], n), getsegments(basisVector[j], m))
                            if !isempty(intersect_nm )
                                if k == 1
                                    push!(interaction_index, CartesianIndex(n, m))
                                end
                                for p in eachindex(basisVector[k].infos)
                                    if !isempty(intersect(intersect_nm , getsegments(basisVector[k], p)))
                                        push!(interaction_index3, CartesianIndex(n,m,p))
                                    end                                    
                                end
                            end
                        end
                    end
                    size_k += length(basisVector[k])
                    block3 = InfoBlock(CartesianIndex(i,j,k), (axes1, axes2, axes3), i == j == k, interaction_index3)
                    blocks3[ib3] = block3
                    ib3+=1
                end
                size_j += length(basisVector[j])
                block = InfoBlock(CartesianIndex(i,j), (axes1, axes2), i == j, interaction_index)
                blocks[ib] = block
                ib+=1
            end
            size_i += length(basisVector[i])
        end
        
        new(basisVector, size, blocks, blocks3, cumul_index)
    end
end

@inline bottom_type(cb::CombineShortPolynomialBasis)= bottom_type(first(cb))

@inline Base.length(cb::CombineShortPolynomialBasis) = cb.size
@inline Base.eachindex(cb::CombineShortPolynomialBasis) = 1:cb.size
@inline Base.first(cb::CombineShortPolynomialBasis) = cb.basisVector[1]
@inline getbasis(cb::CombineShortPolynomialBasis, i::Int) = cb.basisVector[i]
@inline getblocks(cb::CombineShortPolynomialBasis) = cb.blocks
@inline getblocks3(cb::CombineShortPolynomialBasis) = cb.blocks3

function find_basis(cb::CombineShortPolynomialBasis, i::Int)
    @assert i ≤ length(cb)
    ib = 1
    while (ib < length(cb.cumul_index)) && (cb.cumul_index[ib+1] ≤ i) 
        ib += 1
    end
    (ib, i-cb.cumul_index[ib]+1)
end

########################################################################################
#                              Generation of FEM Matrices
########################################################################################

# Mass matrix

function mass_matrix(cb::CombineShortPolynomialBasis)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    fill_mass_matrix!(cb, A)
    A
end

function fill_mass_matrix!(cb::CombineShortPolynomialBasis, A)
    for b ∈ getblocks(cb)
        @views ABlock = A[getrangerow(b), getrangecolumn(b)]
        if isdiagonal(b)
            fill_mass_matrix!(getbasis(cb, _getindex(b,1)), ABlock)
        else
            fill_mass_matrix!(getbasis(cb, _getindex(b,1)), getbasis(cb, _getindex(b,2)), b.interaction_index, ABlock)
            @views ABlockT = A[getrangecolumn(b), getrangerow(b)]
            @. ABlockT = ABlock'
        end
    end
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
        #@inbounds A[I[1], I[2]] *= getnormalization(spb1, I[1]) 
        #@inbounds A[I[1], I[2]] *= getnormalization(spb2, I[2])
    end
end

# Stiffness matrix

function stiffness_matrix(cb::CombineShortPolynomialBasis)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    fill_stiffness_matrix!(cb, A)
    A
end

function fill_stiffness_matrix!(cb::CombineShortPolynomialBasis, A)
    for b ∈ getblocks(cb)
        @views ABlock = A[getrangerow(b), getrangecolumn(b)]
        if isdiagonal(b)
            fill_stiffness_matrix!(getbasis(cb, _getindex(b,1)), ABlock)
        else
            fill_stiffness_matrix!(getbasis(cb, _getindex(b,1)), getbasis(cb, _getindex(b,2)), b.interaction_index, ABlock)
            @views ABlockT = A[getrangecolumn(b), getrangerow(b)]
            @. ABlockT = ABlock'
        end
    end
end

function fill_stiffness_matrix!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, interaction_index::Vector{CartesianIndex{2}}, A)
    for I ∈ interaction_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb1, I[1]), getsegments(spb2, I[2]))
            P = getderivpolynomial(spb1, I[1], i)
            Q = getderivpolynomial(spb2, I[2], j)
            ϕ = getshift(spb1, I[1], i)
            dϕ = ϕ[1]
            @inbounds A[I[1], I[2]] += dϕ * scalar_product(P, Q, spb1.elements.binf, spb1.elements.bsup)
        end
        #@inbounds A[I[1], I[2]] *= getnormalization(spb1, I[1]) 
        #@inbounds A[I[1], I[2]] *= getnormalization(spb2, I[2])
    end
end


# Weight Mass matrix

function weight_mass_matrix(cb::CombineShortPolynomialBasis, weight)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    fill_weight_mass_matrix!(cb, weight, A)
    A
end

function fill_weight_mass_matrix!(cb::CombineShortPolynomialBasis, weight, A)
    for b ∈ getblocks(cb)
        @views ABlock = A[getrangerow(b), getrangecolumn(b)]
        if isdiagonal(b)
            fill_weight_mass_matrix!(getbasis(cb, _getindex(b,1)), weight, ABlock)
        else
            fill_weight_mass_matrix!(getbasis(cb, _getindex(b,1)), getbasis(cb, _getindex(b,2)), weight, b.interaction_index, ABlock)
            @views ABlockT = A[getrangecolumn(b), getrangerow(b)]
            @. ABlockT = ABlock'
        end
    end
end

function weight_mass_matrix(cb::CombineShortPolynomialBasis, n::Int)
    weight_mass_matrix(cb, Monomial(n))
end

function fill_weight_mass_matrix!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, weight, interaction_index::Vector{CartesianIndex{2}}, A)
    for I ∈ interaction_index
        for (i,j) ∈ intersection_with_indices(getsegments(spb1, I[1]), getsegments(spb2, I[2]))
            P = getpolynomial(spb1, I[1], i)
            Q = getpolynomial(spb2, I[2], j)
            invϕ = getinvshift(spb1, I[1], i)
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2]] += dinvϕ * weight_scalar_product(P, Q, weight, spb1.elements.binf, spb1.elements.bsup, invϕ)
        end
        #@inbounds A[I[1], I[2]] *= getnormalization(spb1, I[1]) 
        #@inbounds A[I[1], I[2]] *= getnormalization(spb2, I[2])
    end
end

# Weight Mass vector

function weight_mass_vector(cb::CombineShortPolynomialBasis, weight)
    T = bottom_type(first(cb))
    A = zeros(T, length(cb))
    fill_weight_mass_vector(cb, weight, A)
    A
end

function fill_weight_mass_vector(cb::CombineShortPolynomialBasis, weight, A)
    for b ∈ getblocks(cb)
        if isdiagonal(b)
            @views ABlock = A[getrangerow(b)]
            fill_weight_mass_vector!(getbasis(cb, _getindex(b,1)), weight, ABlock)
        end
    end
end

# Weight Mass tensor

function weight_mass_3tensor(cb::CombineShortPolynomialBasis, weight)
    T = bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb), length(cb)))
    fill_weight_mass_3tensor!(cb, weight, A)
    A
end

function fill_weight_mass_3tensor!(cb::CombineShortPolynomialBasis, weight, A)
    for b ∈ getblocks3(cb)
        @views ABlock = A[getaxes(b,1), getaxes(b,2), getaxes(b,3)]
        if isdiagonal(b)
            fill_weight_mass_3tensor!(getbasis(cb, _getindex(b,1)), weight, ABlock)
        else
            fill_weight_mass_3tensor!(getbasis(cb, _getindex(b,1)), getbasis(cb, _getindex(b,2)), getbasis(cb, _getindex(b,3)), weight, b.interaction_index, ABlock)
            @views ABlockT = A[getaxes(b,1), getaxes(b,3), getaxes(b,2)]
            ABlockT .= permutedims(ABlock, (1,3,2))
            @views ABlockT = A[getaxes(b,2), getaxes(b,1), getaxes(b,3)]
            ABlockT .= permutedims(ABlock, (2,1,3))
            @views ABlockT = A[getaxes(b,2), getaxes(b,3), getaxes(b,1)]
            ABlockT .= permutedims(ABlock, (2,3,1))
            @views ABlockT = A[getaxes(b,3), getaxes(b,1), getaxes(b,2)]
            ABlockT .= permutedims(ABlock, (3,1,2))
            @views ABlockT = A[getaxes(b,3), getaxes(b,2), getaxes(b,1)]
            ABlockT .= permutedims(ABlock, (3,2,1))
        end
    end
end

function fill_weight_mass_3tensor!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, spb3::ShortPolynomialBasis, weight, interaction_index::Vector{CartesianIndex{3}}, A)
    for I ∈ interaction_index
        for (i,j,k) ∈ intersection_with_indices(getsegments(spb1, I[1]), getsegments(spb2, I[2]), getsegments(spb3, I[3]))
            P = getpolynomial(spb1, I[1], i)
            Q = getpolynomial(spb2, I[2], j)
            L = getpolynomial(spb3, I[3], k)
            invϕ = getinvshift(spb1, I[1], i)
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2], I[3]] += dinvϕ * weight_scalar_product(P, Q, L, weight, spb1.elements.binf, spb1.elements.bsup, invϕ)
        end
        #@inbounds A[I[1], I[2], I[3]] *= getnormalization(spb1, I[1]) 
        #@inbounds A[I[1], I[2], I[3]] *= getnormalization(spb2, I[2])
        #@inbounds A[I[1], I[2], I[3]] *= getnormalization(spb3, I[3])
    end
end

function build_basis(cb::CombineShortPolynomialBasis, i::Int)
    (ib, iib) = find_basis(cb, i)
    spb = getbasis(cb, ib)
    build_basis(spb, iib)
end

function eval_basis(cb::CombineShortPolynomialBasis, i::Int, x)
    (ib, iib) = find_basis(cb, i)
    spb = getbasis(cb, ib)
    eval_basis(spb, iib, x)
end

function build_on_basis(cb::CombineShortPolynomialBasis, coeffs)
    @assert eachindex(coeffs) == eachindex(cb) 
    poly = coeffs[firstindex(coeffs)] * build_basis(cb, firstindex(coeffs))
    for i ∈ eachindex(cb)[2:end]
        poly += coeffs[i] * build_basis(cb, i)
    end
    poly
end