using KohnShamResolution
using TimerOutputs
using UnPack

using KohnShamResolution:   init, loopheader!, loopfooter!, makesolution,
                            prepare_eigenvalue_problem!, find_orbital!, aufbau!, density!, update_density!,
                            compute_total_energy, compute_kinetic_energy, compute_coulomb_energy, compute_hartree_energy, isthereExchangeCorrelation, compute_exchangecorrelation_energy 



to = TimerOutput()

# Creation of the model
z = 10
N = 10
KM = KohnShamExtended(z = z, N = N)

# Choice of the method
method = ODA(0.3)

# Discretization 
lₕ = 1
Nₕ = 80
Rmin = 0
Rmax = 80

# INITIALIZATION 

@timeit to "Create mesh" m = linmesh(Rmin, Rmax, Nₕ)

@timeit to "Create basis" basis = P1IntLegendreGenerator(m; ordermax = 5)

@timeit to "init Discretization" discretization = LDADiscretization(lₕ, basis, m)

@timeit to "Init Solver" solver = KohnShamResolution.init(KM, discretization, method; scftol = 1e-3, hartree = false, logconfig = LogConfig(orbitals_energy = true))

# SOLVE FUNCTION

@timeit to "Loop header" loopheader!(solver)


@timeit to "PerformStep" begin

    for i ∈ 1:50
        @unpack model, opts, energies, cache = solver
        @unpack D, Dprev, U, ϵ, n = cache
    
        # STEP 1 : PREPARE THE EIGENVALUE PROBLEM
        @timeit to "Prepare eigenvalue problem" prepare_eigenvalue_problem!(discretization, model, Dprev, opts.hartree)

        # STEP 2 : FIND ORBITALS AND CORRESPONFING ENERGIES
        @timeit to "Find orbital" find_orbital!(discretization, model, U, ϵ)

        # STEP 3 : FILL THE OCCUPATION NUMBER MATRIX ACCORDINGLY WITH THE AUFBAU PRINCIPLE
        @timeit to "aufbau" aufbau!(cache, solver)

        if !cache.flag_degen

            # STEP 4 : COMPUTE A GUESS DENSITY
            @timeit to "density computation" density!(discretization, U, n, D)

            # STEP 5 : COMPUTE ALL ENERGIES
            @timeit to "compute energy" begin 
            energies[:Etot] = compute_total_energy(discretization, model, D, n, ϵ)
            energies[:Ekin] = compute_kinetic_energy(discretization, U, n)
            energies[:Ecou] = compute_coulomb_energy(discretization, U, n)
            energies[:Ehar] = compute_hartree_energy(discretization, D)
            !isthereExchangeCorrelation(model) || (energies[:Eexc] = compute_exchangecorrelation_energy(discretization, model, D))
            end
        end

      # STEP  6 : COMPUTE THE NEW DENSITY
      @timeit to "update density" update_density!(cache, method, solver)  
    end
end

@timeit to "Loop footer" loopfooter!(solver)

@timeit to "Make Solution" makesolution(solver, "")

original_stdout = stdout
output_file = open("tests/Performance/solver/lda_cda.txt", "a")
redirect_stdout(output_file)

print(to)

redirect_stdout(original_stdout)
close(output_file)
println("Performance finished")
