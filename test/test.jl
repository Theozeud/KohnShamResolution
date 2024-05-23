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

# Create a mesh
m = mesh(0.0, x->exp(0.4*x), 10)


# Solve


# Plot Results

