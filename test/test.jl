using KohnShamResolution

# Creation of the model
_z = 2
_N = 10

KM = KohnShamExtended(z = _z,N = _N)

# Choice of the method
method = ODA()

# Discretization 
lₕ = 2
m = mesh(1:5)
basis = HatBasis(m)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solve
@time groundstate(KM, D, method; tol = 0.1)

# Plot Results

