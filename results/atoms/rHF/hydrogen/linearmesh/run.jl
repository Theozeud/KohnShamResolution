include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 3,
                nh              = 5, 
                method          = ODA(0.8), 
                model           = ReducedHartreeFock(21, 21), 
                Rmax            = 300, 
                Nmesh           = 30,
                typemesh        = expmesh, 
                optsmesh        = (s = 1.5,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 15,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-9,
                maxiter         = 100,
                hartree         = true,
                degen_tol       = 1e-3,
                logconfig       = logconfig,
                verbose         = 1)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])