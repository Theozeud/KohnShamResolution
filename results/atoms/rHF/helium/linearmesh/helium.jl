include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# MODEL
z = 2
N = 2
model = ReducedHartreeFock(z, N)

# TYPE
T = Float64

# MESH
Rmin = 0.0
Rmax = 40.0
Nmesh = 40
mesh = linmesh(Rmin, Rmax, Nmesh; T = T)

# BASIS
basis = ShortP1IntLegendreBasis(mesh, T; ordermax = 5)

# DISCRETIZATION
lh = 0
discretization = KohnShamRadialDiscretization(lh, basis, mesh)

# SCF METHOD
method = CDA(0.7)
scftol = 1e-14 
maxiter = 100

# LOG CONFIG
logconfig = LogConfig()

# RESOLUTION
@time sol = groundstate(model, discretization, method; 
        scftol = scftol, 
        maxiter = maxiter, 
        logconfig = logconfig, 
        hartree = true)


plot_stopping_criteria([sol])