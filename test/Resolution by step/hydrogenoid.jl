using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 4

# Discretization 
Rmin = 0
Rmax = 10
Nmesh = 5
m = linmesh(Rmin, Rmax, Nmesh)
basis = ShortP1IntLegendreBasis(m; left = false, right = false, normalize = true, ordermin = 2, ordermax = 2)

# Solve
deriv_basis = deriv(basis)
A   = mass_matrix(deriv_basis)
M₀  = mass_matrix(basis)
M₋₁ = weight_mass_matrix(basis, -1)
M₋₂ = weight_mass_matrix(basis, -2)


H = 1/2 * A  - z .* M₋₁
ϵ, U = eigen(H,M₀)

