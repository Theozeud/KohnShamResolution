using KohnShamResolution

# Creation of the model
_z = 1
_N = 1

KM = KohnShamExtended(z = _z,N = _N)

# Choice of the method
method = ODA()

# Choice of the Solver Options
lₕ = 10
Nₕ = 10

# Solve
groundstate(KM, method; ci = 0.0, lₕ = 2, Nₕ = 2, ε = 0.1)

# Plot Results

