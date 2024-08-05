using KohnShamResolution

# Type of numbers
T = Float64

# Choice of the method
method = ConstantODA(T(0.5))

# Model
z = 1
N = 1
KM = KohnShamExtended(z = z, N = N)

# Discretization 
lₕ = 0

Rmin = 0
Rmax = 200
Nmesh = 80
m = linmesh(Rmin, Rmax, Nmesh)

# Choice of the basis
basis = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 2, normalize = true)

# Final Discretization
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solution
sol = groundstate(KM, D, method; tol = 1e-5, hartree = true, maxiter = 10, potential = :pde)