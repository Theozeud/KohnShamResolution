using KohnShamResolution
using LinearAlgebra
using Plots

# Discretization 
Nmesh = 11
Rmin = 0
Rmax = 10
m = linmesh(Rmin, Rmax, Nmesh)
basis = ShortP1Basis(m; left = false, right = false, normalize = true)

# Solve

deriv_basis = deriv(basis)

@time mass_matrix(deriv_basis)
@time mass_matrix(basis)
@time M₋₁ = weight_mass_matrix(basis, -1)
@time M₋₂ = weight_mass_matrix(basis, -2)
