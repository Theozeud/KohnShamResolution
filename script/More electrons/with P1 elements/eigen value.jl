using KohnShamResolution

# Choice of the method
method = ODA()

# Model
z = 8
N = 4
KM = KohnShamExtended(z = z, N = N)

# Discretization 
lₕ = 4
cutting_pre = 10

Rmin = 0
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, 100; z = 1/z)

# Choice of the basis
basis = HatBasis(m; left = false, right = false)

# Final Discretization
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solution
sol = groundstate(KM, D, method; tol = 1e-4, hartree = false, maxiter = 1000, quad_reltol = 1e-5, quad_abstol = 1e-5)