using KohnShamResolution

# Creation of the model
_z = 1
_N = 1

KM = KohnShamExtended(z = _z,N = _N)

# Choice of the method
method = ODA()

# Choice of the Solver Options
lₕ = 1
Nₕ = 1
basis = LaurentPolynomialBasis([Monomial(1)])

# Solve
groundstate(KM, KohnShamSphericalDiscretization(lₕ,Nₕ,basis), method; tol = 0.1)

# Plot Results

