# UTILS FOR INTERSECTION
function find_intersection_indices(A::Vector{Int}, B::Vector{Int})
    intersect_elements = filter(x -> x != 0, intersect(A, B))
    return [(findfirst(x -> x == el, A), findfirst(x -> x == el, B)) for el in intersect_elements]
end

function find_intersection_indices(A::Vector{Int}, B::Vector{Int}, C::Vector{Int})
    intersect_elements = filter(x -> x != 0, intersect(A, B, C))
    return  [(findfirst(x -> x == el, A), findfirst(x -> x == el, B), findfirst(x -> x == el, C)) for el in intersect_elements]
end

######################       GENERATION OF FEM MATRICES       ########################

# MASS MATRIX
function mass_matrix(pb::PolynomialBasis)
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_mass_matrix!(pb, A)
    A
end

function sparse_mass_matrix(pb::PolynomialBasis)
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_mass_matrix!(pb, A)
    A
end

function fill_mass_matrix!(pb::PolynomialBasis, A::AbstractMatrix{<:Real})
    @unpack generators, mesh = pb
    for I ∈ pb.matrix_fill_indices
        for (i,j) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:])
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            invϕ = getinvshift(pb, I[1], i)
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2]] += dinvϕ * fast_scalar_product(P, Q, pb.generators.binf, pb.generators.bsup)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

# STIFFNESS MATRIX
function stiffness_matrix(pb::PolynomialBasis)
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_stiffness_matrix!(pb, A)
    A
end

function sparse_stiffness_matrix(pb::PolynomialBasis)
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_stiffness_matrix!(pb, A)
    A
end

function fill_stiffness_matrix!(pb::PolynomialBasis, A::AbstractMatrix{<:Real})
    @unpack generators, mesh = pb
    for I ∈ pb.matrix_fill_indices
        for (i,j) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:])
            P = getderivgenerator(pb, I[1], i)
            Q = getderivgenerator(pb, I[2], j)
            ϕ = getshift(pb, I[1], i)
            dϕ = ϕ[1]
            @inbounds A[I[1], I[2]] += dϕ * fast_scalar_product(P, Q, pb.generators.binf, pb.generators.bsup)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
end

# WEIGHT MASS MATRIX
function weight_mass_matrix(pb::PolynomialBasis, weight)
    T = eltype(pb)
    A = zeros(T, pb.size, pb.size)
    fill_weight_mass_matrix!(pb, weight, A)
    A
end

function sparse_weight_mass_matrix(pb::PolynomialBasis, weight)
    T = eltype(pb)
    A = spzeros(T, pb.size, pb.size)
    fill_weight_mass_matrix!(pb, weight, A)
    A
end

function weight_mass_matrix(pb::PolynomialBasis, n::Int)
    weight_mass_matrix(pb, Monomial(n))
end

function fill_weight_mass_matrix!(pb::PolynomialBasis, weight, A::AbstractMatrix{<:Real})
    for I ∈ pb.matrix_fill_indices
        for (i,j) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:])
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            invϕ = getinvshift(pb, I[1], i)
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2]] += dinvϕ * weight_scalar_product(P, Q, weight, pb.generators.binf, pb.generators.bsup, invϕ)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

function fill_weight_mass_matrix!(pb::PolynomialBasis, n::Int, A::AbstractMatrix{<:Real})
    fill_weight_mass_matrix!(pb::PolynomialBasis, Monomial(n), A)
end

# WEIGHT MASS VECTOR
function weight_mass_vector(pb::PolynomialBasis, weight)
    T = eltype(pb)
    A = zeros(T, pb.size)
    fill_weight_mass_vector!(pb, weight, A)
    A
end

function sparse_weight_mass_vector(pb::PolynomialBasis, weight)
    T = eltype(pb)
    A = spzeros(T, pb.size)
    fill_weight_mass_vector!(pb, weight, A)
    A
end

function fill_weight_mass_vector!(pb::PolynomialBasis, weight, A::AbstractMatrix{<:Real})
    for i ∈ eachindex(pb)
        for j ∈ axes(pb.indices_generators,2)
            P = getgenerator(pb, i, j)
            invϕ = getinvshift(pb, i, j)
            dinvϕ = invϕ[1]
            @inbounds A[i] += dinvϕ * weight_scalar_product(P, weight, pb.generators.binf, pb.generators.bsup, invϕ)
        end
    end
end

# WEIGHT MASS TENSOR
function weight_mass_3tensor(pb::PolynomialBasis, weight)
    T = eltype(pb)
    A = zeros(T, pb.size, pb.size, pb.size)
    fill_weight_mass_3tensor!(pb, weight, A)
    A
end

function weight_mass_3tensor(pb::PolynomialBasis, n::Int)
    weight_mass_3tensor(pb, Monomial(n))
end


function fill_weight_mass_3tensor!(pb::PolynomialBasis, weight, 
                                   A::AbstractArray{<:Real})
    for I ∈ pb.tensor_fill_indices
        for (i,j,k) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:], pb.indices_cells[I[3],:])
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            L = getgenerator(pb, I[3], k)
            invϕ = getinvshift(pb, I[1], i)
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2], I[3]] += dinvϕ * weight_scalar_product(P, Q, L, weight, pb.generators.binf, pb.generators.bsup, invϕ)
        end
        @inbounds A[I[3], I[1], I[2]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[2], I[3], I[1]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[2], I[1], I[3]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[3], I[2], I[1]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[1], I[3], I[2]]  = A[I[1], I[2], I[3]]
    end
    nothing
end

function fill_weight_mass_3tensor!(pb::PolynomialBasis, n::Int, A::AbstractArray{<:Real})
    fill_weight_mass_3tensor!(pb, Monomial(n),A)
end


function fill_weight_mass_3tensor!(pb::PolynomialBasis, weight, 
                                   A::Dict{Tuple{Int64, Int64, Int64}, <:Real})
    for I ∈ pb.tensor_fill_indices
        for (i,j,k) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:], pb.indices_cells[I[3],:])
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            L = getgenerator(pb, I[3], k)
            invϕ = getinvshift(pb, I[1], i)
            dinvϕ = invϕ[1]
            if haskey(A, (I[1], I[2], I[3]))
                A[(I[1], I[2], I[3])] += dinvϕ * weight_scalar_product(P, Q, L, weight, pb.generators.binf, pb.generators.bsup, invϕ)
            else
                A[(I[1], I[2], I[3])]  = dinvϕ * weight_scalar_product(P, Q, L, weight, pb.generators.binf, pb.generators.bsup, invϕ)
            end
        end
        A[(I[3], I[1], I[2])]  = A[(I[1], I[2], I[3])]
        A[(I[2], I[3], I[1])]  = A[(I[1], I[2], I[3])]
        A[(I[2], I[1], I[3])]  = A[(I[1], I[2], I[3])]
        A[(I[3], I[2], I[1])]  = A[(I[1], I[2], I[3])]
        A[(I[1], I[3], I[2])]  = A[(I[1], I[2], I[3])]
    end
    nothing
end