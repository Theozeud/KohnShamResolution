using KohnShamResolution
using LinearAlgebra
using TimerOutputs
using Plots

to = TimerOutput()

# Parameters
model   = LSDA(1, 1)
method  = CDA(0.7)
lₕ = 0
Nₕ = 40
Rmin = 0
Rmax = 40
mesh = linmesh(Rmin, Rmax, Nₕ)

# Initialisation

@timeit to "Initialisation" begin

    @timeit to "Create Basis" basis = P1IntLegendreGenerator(mesh; ordermax = 3)

    @timeit to "Create Discretization" discretization = LSDADiscretization(lₕ, basis, mesh)

    @timeit to "Init Solver" solver = KohnShamResolution.init(model, discretization, method; scftol = 1e-3, hartree = false)
end


@timeit to "PerformStep" begin
    @timeit to "Find orbitals" KohnShamResolution.find_orbital!(solver.discretization, solver)

    @timeit to "Aufbau" KohnShamResolution.aufbau!(solver)

    @timeit to "Normalization" KohnShamResolution.normalization!(solver.discretization)

    @timeit to "Compute Density matrix" KohnShamResolution.density_matrix!(solver.discretization)

    @timeit to "SCF Method" KohnShamResolution.update_density!(solver.method, solver)
end

@timeit to "Compute Energy" begin
    @timeit to "Kinetic Energy" KohnShamResolution.compute_kinetic_energy!(discretization,solver)
    @timeit to "Coulomb Energy" KohnShamResolution.compute_coulomb_energy!(discretization,solver)
    @timeit to "Hartree Energy" KohnShamResolution.compute_hartree_energy!(discretization,solver)
    @timeit to "Exch Energy"    KohnShamResolution.compute_exchangecorrelation_energy!(discretization,solver)
    @timeit to "KinCorr Energy" KohnShamResolution.compute_kinetic_correlation_energy!(discretization,solver)
    @timeit to "Total Energy"   KohnShamResolution.compute_total_energy!(discretization,solver)
end

@timeit to "Loop footer" KohnShamResolution.loopfooter!(solver)

@timeit to "Make Solution" KohnShamResolution.makesolution(solver, "")


to