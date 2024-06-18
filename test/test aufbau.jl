using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 4
N = 14

KM = KohnShamExtended(z = z,N = N)
#KM = SlaterXα(z, N)

# Choice of the method
method = ODA()

# Discretization 
lₕ = 2
Rmin = 0
cutting_pre = 10
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, 10)
basis = P1Basis(m; left = false, right = false)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solve
#@time sol = groundstate(KM, D, method; tol = 1e-10, hartree = false)

deriv_basis = deriv(basis)
 
A   = mass_matrix(deriv_basis, Rmin, Rmax)
M₀  = mass_matrix(basis, Rmin, Rmax)
M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

U = zeros(lₕ+1, 8, 8)
ϵ = zeros(lₕ+1, 8)

for l ∈ 0:lₕ
    H = 1/2 * (A + l*(l+1)*M₋₂) - z .* M₋₁
    _ϵ, _U = eigen(H,M₀)
    ϵ[l+1,:] .= _ϵ
    U[l+1,:,:] .= _U
end

display(ϵ)

n = zeros(lₕ+1, 8)
function aufbau!(n, ϵ, N; tol = eps(eltype(ϵ)))
    _l,_n = size(ϵ)
    ϵ_vect = vec(ϵ)
    @show index_sort = sortperm(ϵ_vect)
    degen_matrix = reduce(hcat, [[2*l + 1 for l ∈ 0:_l-1] for i ∈ 1:_n])
    remain = N
    idx = 1
    
    while remain > 0 && idx < length(ϵ) + 1

        A = Int[]  #Stock all index corresponding to a degenerancy
        ϵ_cur = ϵ_vect[index_sort[idx]]
        push!(A, index_sort[idx])
        idx += 1
        while abs(ϵ[index_sort[idx]] - ϵ_cur) < tol
            push!(A, index_sort[idx])
            idx += 1
        end
        @show A
        # Count total degeneracy
        degen = sum(degen_matrix[i] for i in A)
        
        # See what to do depending on the case
        if remain - degen ≥ 0
            for i in A
                n[i] = 2 * degen_matrix[i]
            end
            remain -= degen
        else
            if length(A) == 1
                # First case, if no degeneracy
                n[first(A)] = 2 * remain
            else
                # Second case, if degeneracy
                @error "There is accidental degeneracy but no implementation for this case for the moment."
            end
            break
        end
    end
end

aufbau!(n, ϵ, N)
display(n)