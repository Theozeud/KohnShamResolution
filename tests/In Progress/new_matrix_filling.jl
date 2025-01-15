using KohnShamResolution
using KohnShamResolution: getgenerator, getinvshift, find_intersection_indices, fast_scalar_product
using UnPack
using Base.Threads

mesh = linmesh(0,30,300)
pb = P1IntLegendreGenerator(mesh, Float64; ordermax = 4)

function mass_matrix2(pb::PolynomialBasis, matrix_generator)
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_mass_matrix2!(pb, A, matrix_generator)
    A
end

function matrix_gen(pb)
    matrix_generator = zeros(eltype(pb), length(pb.generators), length(pb.generators))
    for i ∈ eachindex(pb.generators)
        for j ∈ i:length(pb.generators)
            P = pb.generators[i]
            Q = pb.generators[j]
            matrix_generator[i,j] = scalar_product(P, Q, pb.generators.binf, pb.generators.bsup)
            matrix_generator[j,i] = matrix_generator[i,j]
        end
    end
    matrix_generator
end

function fill_mass_matrix2!(pb::PolynomialBasis, A, matrix_generator)
    @unpack generators, mesh = pb
    
    @threads for I ∈ pb.matrix_fill_indices
        for (i,j) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:])
            k1 = pb.indices_generators[I[1],i]
            k2 = pb.indices_generators[I[2],j]
            invϕ = getinvshift(pb, I[1], i)
            dinvϕ = invϕ[1]
            @inbounds A[I[1], I[2]] += dinvϕ * matrix_generator[k1,k2]
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

Mgen = matrix_gen(pb)

@time Mold = mass_matrix(pb)
@time Mnew = mass_matrix2(pb, Mgen)

using Test
@test Mold == Mnew