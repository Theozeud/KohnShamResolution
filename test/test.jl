using KohnShamResolution

# Creation of the model
_z = 1
_N = 1

KM = KohnShamExtended(z = _z,N = _N)

# Choice of the method
method = ODA()

# Discretization 
lₕ = 2
m = mesh(1:10)
basis = HatBasis(m)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solve
@time groundstate(KM, D, method; tol = 0.1)

# Plot Results

