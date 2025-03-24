# UTILS FOR INTERSECTION

function find_intersection_indices(A::Vector{Int}, B::Vector{Int})
    intersect_elements = filter(x -> x != 0, intersect(A, B))
    return [(findfirst(x -> x == el, A), findfirst(x -> x == el, B)) for el in intersect_elements]
end

function find_intersection_indices(A::Vector{Int}, B::Vector{Int}, C::Vector{Int})
    intersect_elements = filter(x -> x != 0, intersect(A, B, C))
    return  [(findfirst(x -> x == el, A), findfirst(x -> x == el, B), findfirst(x -> x == el, C)) for el in intersect_elements]
end

#####################################################################
#                          OVERLAP MATRIX
#####################################################################

function mass_matrix(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = ExactIntegration())
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_mass_matrix!(pb, A; weight = weight, method = method)
    A
end

function sparse_mass_matrix(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = ExactIntegration())
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_mass_matrix!(pb, A; weight = weight, method = method)
    A
end

function mass_matrix(pb::PolynomialBasis, n::Int; method::IntegrationMethod = ExactIntegration())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_matrix(pb; weight = InvX(), method = method)
    elseif n == -2
        mass_matrix(pb; weight = InvX2(), method = method)
    end
end

function sparse_mass_matrix(pb::PolynomialBasis, n::Int; method::IntegrationMethod = ExactIntegration())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        sparse_mass_matrix(pb; weight = InvX(), method = method)
    elseif n == -2
        sparse_mass_matrix(pb; weight = InvX2(), method = method)
    end
end


function fill_mass_matrix!( pb::PolynomialBasis, 
                            A::AbstractMatrix{<:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = ExactIntegration())
    @unpack generators, mesh = pb
    for I ∈ pb.matrix_fill_indices
        for (i,j) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:])
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            ϕ = getshift(pb, I[1], i)
            invϕ = getinvshift(pb, I[1], i)
            a = mesh[pb.indices_cells[I[1],i]]
            b = mesh[pb.indices_cells[I[1],i]+1]
            intdata = IntegrationData(weight,
                                      (P,Q),
                                      ϕ,
                                      invϕ,
                                      a,
                                      b,
                                      pb.generators.binf,
                                      pb.generators.bsup,
                                      method)
            @inbounds A[I[1], I[2]] += swsp(intdata)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end


function fill_mass_matrix!(pb::PolynomialBasis, n::Int, A::AbstractArray{<:Real}; method::IntegrationMethod = ExactIntegration())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        fill_mass_matrix!(pb, A; weight = InvX(), method = method)
    elseif n == -2
        fill_mass_matrix!(pb, A; weight = InvX2(), method = method)
    end
    nothing
end

#####################################################################
#                          STIFFNESS MATRIX
#####################################################################

function stiffness_matrix(pb::PolynomialBasis; method::IntegrationMethod = ExactIntegration())
    @unpack size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_stiffness_matrix!(pb, A; method = method)
    A
end

function sparse_stiffness_matrix(pb::PolynomialBasis; method::IntegrationMethod = ExactIntegration())
    @unpack size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_stiffness_matrix!(pb, A; method = method)
    A
end

function fill_stiffness_matrix!(pb::PolynomialBasis, A::AbstractMatrix{<:Real}; method::IntegrationMethod = ExactIntegration())
    @unpack mesh = pb
    for I ∈ pb.matrix_fill_indices
        for (i,j) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:])
            P = getderivgenerator(pb, I[1], i)
            Q = getderivgenerator(pb, I[2], j)
            ϕ = getshift(pb, I[1], i)
            invϕ = getinvshift(pb, I[1], i)
            a = mesh[pb.indices_cells[I[1],i]]
            b = mesh[pb.indices_cells[I[1],i]+1]
            intdata = IntegrationData(NoWeight(),
                                      (P,Q),
                                      ϕ,
                                      invϕ,
                                      a,
                                      b,
                                      pb.generators.binf,
                                      pb.generators.bsup,
                                      method)
            @inbounds A[I[1], I[2]] += ϕ[1] * swsp(intdata)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
end

#####################################################################
#                             MASS VECTOR
#####################################################################

function mass_vector(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = ExactIntegration())
    T = eltype(pb)
    A = zeros(T, pb.size)
    fill_mass_vector!(pb, A; weight = weight, method = method)
    A
end

function sparse_mass_vector(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = ExactIntegration())
    T = eltype(pb)
    A = spzeros(T, pb.size)
    fill_mass_vector!(pb, A; weight = weight, method = method)
    A
end

function fill_mass_vector!(pb::PolynomialBasis, A::AbstractMatrix{<:Real}; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = ExactIntegration())
    @unpack mesh = pb
    for i ∈ eachindex(pb)
        for j ∈ axes(pb.indices_generators,2)
            P = getgenerator(pb, I[1], i)
            ϕ = getshift(pb, I[1], i)
            invϕ = getinvshift(pb, I[1], i)
            a = mesh[pb.indices_cells[I[1],i]]
            b = mesh[pb.indices_cells[I[1],i]+1]
            intdata = IntegrationData(weight,
                                      (P,),
                                      ϕ,
                                      invϕ,
                                      a,
                                      b,
                                      pb.generators.binf,
                                      pb.generators.bsup,
                                      method = method)
            @inbounds A[I[1], I[2]] += swsp(intdata)
        end
    end
end

#####################################################################
#                          MASS TENSOR
#####################################################################

function mass_tensor(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = ExactIntegration())
    T = eltype(pb)
    A = zeros(T, pb.size, pb.size, pb.size)
    fill_mass_tensor!(pb, A; weight = weight, method = method)
    A
end

function mass_tensor(pb::PolynomialBasis, n::Int; method::IntegrationMethod = ExactIntegration())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_tensor(pb; weight = InvX(), method = method)
    elseif n == -2
        mass_tensor(pb; weight = InvX2(), method = method)
    end
end

function fill_mass_tensor!( pb::PolynomialBasis, 
                            A::AbstractArray{<:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = ExactIntegration())
    @unpack mesh = pb
    for I ∈ pb.tensor_fill_indices
        for (i,j,k) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:], pb.indices_cells[I[3],:])
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            L = getgenerator(pb, I[3], k)
            ϕ = getshift(pb, I[1], i)
            invϕ = getinvshift(pb, I[1], i)
            a = mesh[pb.indices_cells[I[1],i]]
            b = mesh[pb.indices_cells[I[1],i]+1]
            intdata = IntegrationData(weight,
                                      (P,Q,L),
                                      ϕ,
                                      invϕ,
                                      a,
                                      b,
                                      pb.generators.binf,
                                      pb.generators.bsup,
                                      method)
            @inbounds A[I[1], I[2]] += swsp(intdata)
        end
        @inbounds A[I[3], I[1], I[2]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[2], I[3], I[1]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[2], I[1], I[3]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[3], I[2], I[1]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[1], I[3], I[2]]  = A[I[1], I[2], I[3]]
    end
    nothing
end

function fill_mass_tensor!(pb::PolynomialBasis, 
                           n::Int, 
                           A::Union{AbstractArray{<:Real},Dict{Tuple{Int64, Int64, Int64}, <:Real}};
                           method::IntegrationMethod = ExactIntegration())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        fill_mass_tensor!(pb, A; weight = InvX(), method = method)
    elseif n == -2
        fill_mass_tensor!(pb, A; weight = InvX2(), method = method)
    end
    nothing
end

function fill_mass_tensor!( pb::PolynomialBasis,
                            A::Dict{Tuple{Int64, Int64, Int64}, <:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = ExactIntegration())
    @unpack mesh = pb
    for I ∈ pb.tensor_fill_indices
        for (i,j,k) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:], pb.indices_cells[I[3],:])
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            L = getgenerator(pb, I[3], k)
            ϕ = getshift(pb, I[1], i)
            invϕ = getinvshift(pb, I[1], i)
            a = mesh[pb.indices_cells[I[1],i]]
            b = mesh[pb.indices_cells[I[1],i]+1]
            intdata = IntegrationData(weight,
                                      (P,Q,L),
                                      ϕ,
                                      invϕ,
                                      a,
                                      b,
                                      pb.generators.binf,
                                      pb.generators.bsup,
                                      method)
            if haskey(A, (I[1], I[2], I[3]))
                A[(I[1], I[2], I[3])] += swsp(intdata)
            else
                A[(I[1], I[2], I[3])]  = swsp(intdata)
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