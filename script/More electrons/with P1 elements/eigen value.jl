using KohnShamResolution

# Choice of the method
method = ODA()

# Model
z = 2
N = 2
KM = KohnShamExtended(z = z, N = N)

# Discretization 
lₕ = 0
cutting_pre = 10

Rmin = 0
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, 120; z = 1/z)

# Choice of the basis
basis = HatBasis(m; left = false, right = false)

# Final Discretization
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solution
sol = groundstate(KM, D, method; tol = 1e-4, hartree = true, maxiter = 1000, quad_reltol = 1e-5, quad_abstol = 1e-5)