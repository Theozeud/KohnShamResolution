using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 32
N = 32

# Choice of the method
method = CDA(1)

# Discretization 
lₕ = 3
Rmin = 0
Rmax =  30
Nmesh = 200
m = linmesh(Rmin, Rmax, Nmesh)
basis = P1IntLegendreGenerator(m; ordermax = 3)
D = LDADiscretization(lₕ, basis, m)

# Solve
A   = stiffness_matrix(basis)
M₀  = mass_matrix(basis)
M₋₁ = weight_mass_matrix(basis, -1)
M₋₂ = weight_mass_matrix(basis, -2)

U = zeros(lₕ+1, length(basis), length(basis))
ϵ = zeros(lₕ+1, length(basis))

for l ∈ 0:lₕ
    H = 1/2 * (A + l*(l+1)*M₋₂) - z .* M₋₁
    _ϵ, _U = eigen(H,M₀)
    ϵ[l+1,:] .= _ϵ
    U[l+1,:,:] .= _U
end

display(ϵ)

n = zeros(lₕ+1, length(basis))

function aufbau!(n, ϵ, N; tol = eps(eltype(ϵ)))

    index_sort = sortperm(vec(ϵ))
    
    function degeneracy(idx)
        l = rem(idx - 1, size(ϵ, 1)) 
        return 4 * l + 2
    end

    remain = N
    idx = 1
    while remain > 0 && idx ≤ length(index_sort)

        indices_degen = Int[]  
        push!(indices_degen, index_sort[idx])
        idx += 1
        while idx ≤ length(index_sort) && abs(ϵ[index_sort[idx]] - ϵ[first(indices_degen)]) < tol
            push!(indices_degen, index_sort[idx])
            idx += 1
        end

        # Count total degeneracy
        total_degen = sum(degeneracy(i) for i in indices_degen)

        if remain - total_degen≥ 0
            for i in indices_degen
                n[i] = degeneracy(i)
            end
            remain -= total_degen
        else
            if length(indices_degen) == 1
                # First case, if no degeneracy
                n[first(indices_degen)] = remain
            else
                # Second case, if degeneracy
                @error "There is accidental degeneracy but no implementation for this case for the moment."
            end
            remain = zero(remain)
        end
    end
end

@time aufbau!(n, ϵ, N)
display(n)