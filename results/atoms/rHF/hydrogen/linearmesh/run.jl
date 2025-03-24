include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0,
                nh              = 2, 
                method          = ODA(0.8), 
                model           = ReducedHartreeFock(2, 2), 
                Rmax            = 100, 
                Nmesh           = 40,
                typemesh        = expmesh, 
                optsmesh        = (s = 2.0,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 10,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-9,
                maxiter         = 40,
                hartree         = true,
                degen_tol       = 1e-10,
                logconfig       = logconfig,
                verbose         = 3)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])