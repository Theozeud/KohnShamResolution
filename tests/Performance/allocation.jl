using KohnShamResolution
using LinearAlgebra
using TimerOutputs

using KohnShamResolution: init, find_orbital!, aufbau!, normalization!, 
density_matrix!, compute_energy!, update_density!, update_solver!, loopfooter!, 
makesolution, performstep!,  loopheader!

to = TimerOutput()

# Creation of the model
z = 1
N = 1
KM = KohnShamExtended(z = z, N = N)

# Choice of the method
M = CDA(0.0)

# Discretization 
lₕ = 0
Nₕ = 100
Rmin = 0
Rmax = 60

@timeit to "Create mesh" m = linmesh(Rmin, Rmax, Nₕ)

@timeit to "Create basis" basis = ShortP1IntLegendreBasis(m; left = false, right = false, ordermin = 2, ordermax = 3)

@timeit to "init Discretization" D = LDADiscretization(lₕ, basis, m)

@timeit to "Init Solver" solver = KohnShamResolution.init(KM, D, M; scftol = 1e-3, hartree = false)



@timeit to "Loop header" loopheader!(solver)

@timeit to "PerformStep" begin

    for i ∈ 1:10
    # STEP 1 : Resolution of the generalized eigenvalue problem to find atomic orbitals and corresonding energies
    @timeit to "Find orbital" find_orbital!(solver.discretization, solver)

    # STEP 2 : Build the n matrix using the Aufbau principle
    @timeit to "aufbau" aufbau!(solver)
    
    # STEP 3 : Normaization of eigenvectors
    @timeit to "normalization" normalization!(solver.discretization, solver)
    
    # STEP 4 : Compute density star
    @timeit to "density matrix computation" density_matrix!(solver.discretization, solver)
    
    # STEP 5 : Compute energy
    @timeit to "compute energy" compute_energy!(solver.discretization, solver)
    
    # STEP 6 : Compute new density
    @timeit to "update density" update_density!(solver.method, solver)
    
    # STEP 7 : Update Solver
    @timeit to "update solver" update_solver!(solver)
    end
end

@timeit to "Loop footer" loopfooter!(solver)

@timeit to "Make Solution" makesolution(solver, "")


to


