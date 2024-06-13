using KohnShamResolution

# Creation of the model
_z = 0.5
_N = 1

KM = KohnShamExtended(z = _z,N = _N)
#KM = SlaterXα(_z, _N)

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0

m = logmesh(0,50,100)
basis = HatBasis(m; left = false, right = false)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solve
@time sol = groundstate(KM, D, method; tol = 1e-20, hartree = false)