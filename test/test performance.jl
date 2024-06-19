using KohnShamResolution
using LinearAlgebra
using TimerOutputs
using Plots

to = TimerOutput()

# Creation of the model
z = 4
N = 1
KM = KohnShamExtended(z = z, N = N)

# Choice of the method
M = ODA()

# Discretization 
lₕ = 0
Nₕ = 20
Rmin = 0
Rmax = 100

@timeit to "Create mesh" m = logmesh(Rmin, Rmax, Nₕ)
@timeit to "Create basis" basis = BubbleBasis(m; order = 3, left = false, right = false)

D = KohnShamSphericalDiscretization(lₕ, basis, m)

@timeit to "Init Solver" solver = KohnShamResolution.init(KM, D, M; tol = 1e-3, hartree = false)

@timeit to "PerformStep" KohnShamResolution.performstep!(M, solver)

@timeit to "Loop footer" KohnShamResolution.loopfooter!(solver, M)

@timeit to "Make Solution" KohnShamResolution.makesolution(solver)

@show to