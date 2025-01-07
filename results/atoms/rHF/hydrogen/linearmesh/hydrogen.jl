include("../../../../../benchmarktools/atoms/includes.jl")
using KohnShamResolution

# MODEL
z = 1
N = 1
model = ReducedHartreeFock(z, N)

# TYPE
T = Float64

# MESH
Rmin = 0.0
Rmax = 20.0
Nmesh = 300
mesh = linmesh(Rmin, Rmax, Nmesh; T = T)

# BASIS
basis = ShortIntLegendreBasis(mesh, T; normalize = false, ordermax = 2)

# DISCRETIZATION
lh = 0
discretization = KohnShamRadialDiscretization(lh, basis, mesh)

# SCF METHOD
method = CDA(0.8)
scftol = 1e-14 
maxiter = 100

# LOG CONFIG
logconfig = LogConfig()

# RESOLUTION
sol = groundstate(model, discretization, method; scftol = scftol, maxiter = maxiter, logconfig = logconfig)
