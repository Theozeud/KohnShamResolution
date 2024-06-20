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
Nₕ = 200
lₕ = 2
Rmin = 0
cutting_pre = 10
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, Nₕ)
basis = P2Basis(m; left = false, right = false)

# Solve

deriv_basis = deriv(basis)
 
@time A   = mass_matrix(deriv_basis, Rmin, Rmax)
@time M₀  = mass_matrix(basis, Rmin, Rmax)
@time M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
@time M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)
