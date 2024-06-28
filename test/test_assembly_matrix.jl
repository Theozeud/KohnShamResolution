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
basis = P1Basis(m; left = false, right = false)

# Solve

#deriv_basis = deriv(basis)

@time mass_matrix(deriv_basis)
@time mass_matrix(basis)
@time M₋₁ = weight_mass_matrix(basis, -1)
@time M₋₂ = weight_mass_matrix(basis, -2)