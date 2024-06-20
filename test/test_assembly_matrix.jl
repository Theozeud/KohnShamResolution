using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 10
N = 5

KM = KohnShamExtended(z = z,N = N)

# Choice of the method
method = ODA()

# Discretization 
Nₕ = 2000
lₕ = 2
Rmin = 0
cutting_pre = 10
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, Nₕ)
basis = P2Basis(m; left = false, right = false)

# Solve

deriv_basis = deriv(basis)

@time mass_matrix(deriv_basis, Rmin, Rmax)
@time mass_matrix(basis, Rmin, Rmax)
@time M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
@time M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

#=
# test New mass matrix
function generate_cross_index(basis)
    cross_index = CartesianIndex[]
    for i in eachindex(basis)
        for j in 1:i
            if !isempty(intersect(basis[i].index, basis[j].index))
                push!(cross_index, CartesianIndex(i, j))
            end
        end
    end
    return cross_index
end
cross_index = generate_cross_index(basis)

function mass_matrix2(lpb::LaurentPolynomialBasis, a::Real, b::Real)
    T = eltype(first(lpb))
    A = zeros(T, (length(lpb), length(lpb)))
    for I ∈ cross_index
        @inbounds A[I[1],I[2]] = scalar_product(lpb[I[1]], lpb[I[2]], a, b)
        @inbounds A[I[2],I[1]] = A[I[1],I[2]] 
    end
    A
end

@time mass_matrix2(basis, Rmin, Rmax)
nothing
=#