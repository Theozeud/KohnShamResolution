using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 3
N = 3

KM = KohnShamExtended(z = z,N = N)
#KM = SlaterXα(z, N)

# Choice of the method
method = ConstantODA(1)

# Discretization 
lₕ = 1
Rmin = 0
Rmax =  15
Nmesh = 500
m = linmesh(Rmin, Rmax, Nmesh)
basis = ShortP1Basis(m; left = false, right = false, normalize = true)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solve

deriv_basis = deriv(basis)
A   = mass_matrix(deriv_basis)
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

#=
1s -> 2s -> 2p || 3s
function next_allowed_occup(last, resi)
    if last == (1,1)
        return (1,2)
    else
        (l,n) = last
        return [(l+1,n-1),

    push!(next_allowed_occup(resi), 
end
(1,1) -> (1,2)

=#

function aufbau!(n, ϵ, N; tol = eps(eltype(ϵ)))
    ϵ_copy = copy(ϵ)
    _l,_n = size(ϵ)
    ϵ_vect = vec(ϵ_copy)
    @show index_sort = sortperm(ϵ_vect)
    @show degen_matrix = reduce(hcat, [[2*l + 1 for l ∈ 0:_l-1] for i ∈ 1:_n])
    remain = N
    idx = 1
    
    while remain > 0 && idx < length(ϵ) + 1

        A = Int[]  # Stock all index corresponding to a degenerancy
        ϵ_cur = ϵ_vect[index_sort[idx]]
        push!(A, index_sort[idx])
        idx += 1
        while abs(ϵ[index_sort[idx]] - ϵ_cur) < tol
            push!(A, index_sort[idx])
            idx += 1
        end
        @show A
        # Count total degeneracy
        @show degen = sum(degen_matrix[i] for i in A)
        
        # See what to do depending on the case
        @show remain
        @show degen
        if remain - 2*degen ≥ 0
            for i in A
                n[i] = 2*degen_matrix[i]
            end
            remain -= 2*degen
        else
            if length(A) == 1
                # First case, if no degeneracy
                n[first(A)] = remain
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


index = findall(x->x ≠ 0, n)
ϵ_full = ϵ[index]
index_sort = sortperm(ϵ_full)
new_index = index[index_sort]

T = eltype(ϵ)
result = Tuple{String,T,T}[]

function convert_into_l(l::Int)
    if l == 0
        return "s"
    elseif l == 1
        return "p"
    elseif l == 2
        return "d"
    elseif l == 3
        return "f" 
    else
        @error "Not implemented Car yet for l ="*string(l)
    end
end 

for i ∈ new_index
    _n = i[2]+i[1]-1
    l = convert_into_l(i[1]-1)
    push!(result,(string(_n,l),ϵ[i], n[i]))
end

result