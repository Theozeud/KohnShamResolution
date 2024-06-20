using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Float64
Rmin = 0
cutting_pre = 50
Rmax = (1.5  + cutting_pre*log(10))
Nmesh = 100
m = linmesh(Rmin, Rmax, Nmesh)

order = 1
basis = IntLegendreBasis(m, T; order = order, left = false, right = false)

X = LinRange(Rmin, Rmax,10*Nmesh)
plt = plot()
for p ∈ basis
    plot!(X,p.(X))
end
plt

deriv_basis = deriv(basis)
A   = mass_matrix(deriv_basis, Rmin, Rmax)
M₀  = mass_matrix(basis, Rmin, Rmax)
M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

z = 1
l = 0
H = 1/2 * (A + l*(l+1)*M₋₂) - z .* M₋₁

using LinearAlgebra
values, vectors = eigen(H, M₀)

#=
for i∈eachindex(values)
    if imag(values[i]) ≠ 0
        println(i)
    end
end
=#

function check_symmetry(M)
    for I ∈ CartesianIndices(M)
        if M[I[1],I[2]] ≠ M[I[2],I[1]]
            print(I)
        end
    end
end

check_symmetry(M₀)
