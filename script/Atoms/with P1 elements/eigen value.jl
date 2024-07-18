using KohnShamResolution

# Choice of the method
method = ConstantODA(0.8)

# Model
z = 1
N = 1
KM = KohnShamExtended(z = z, N = N)

# Discretization 
lₕ = 0
cutting_pre = 10

Rmin = 0
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = linmesh(Rmin, Rmax, 100)

# Choice of the basis
basis = ShortP1IntLegendreBasis(m; left = false, right = false, ordermin = 2, ordermax = 2)

# Final Discretization
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solution
sol = groundstate(KM, D, method; tol = 1e-4, hartree = false, maxiter = 1000, quad_reltol = 1e-5, quad_abstol = 1e-5)