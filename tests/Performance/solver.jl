using KohnShamResolution
using LinearAlgebra
using TimerOutputs
using Plots

to = TimerOutput()

# Creation of the model
z = 1
N = 1
KM = KohnShamExtended(z = z, N = N)

# Choice of the method
M = CDA(0.0)

# Discretization 
lₕ = 3
Nₕ = 300
Rmin = 0
Rmax = 41

@timeit to "Create mesh" m = linmesh(Rmin, Rmax, Nₕ)

@timeit to "Create basis" basis = ShortP1IntLegendreBasis(m; left = false, right = false, normalize = true, ordermin = 2, ordermax = 3)

@timeit to "init Discretization" D = KohnShamRadialDiscretization(lₕ, basis, m)

@timeit to "Init Solver" solver = KohnShamResolution.init(KM, D, M; scftol = 1e-3, hartree = false)

#@timeit to "PerformStep" KohnShamResolution.performstep!(solver)

@timeit to "PerformStep" begin
    @timeit to "Find orbital" KohnShamResolution.find_orbital!(solver.discretization, solver)

    @timeit to "aufbau" KohnShamResolution.aufbau!(solver)

    @timeit to "density matrix computation" KohnShamResolution.density_matrix!(solver.discretization)

    @timeit to "update density" KohnShamResolution.update_density!(solver.method, solver)

    @timeit to "reset cache" KohnShamResolution. update_solver!(solver)
end

@timeit to "Loop footer" KohnShamResolution.loopfooter!(solver)

@timeit to "Make Solution" KohnShamResolution.makesolution(solver)

#@timeit to "Build Density in Solution" KohnShamResolution.build_density2!(solver.discretization, solver.D)

to