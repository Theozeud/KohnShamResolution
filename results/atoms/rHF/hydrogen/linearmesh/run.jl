include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0,
                nh              = 2, 
                method          = ODA(0.8), 
                model           = ReducedHartreeFock(1, 1), 
                Rmax            = 80, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 3,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-12,
                maxiter         = 40,
                hartree         = true,
                degen_tol       = 1e-10,
                logconfig       = logconfig,
                verbose         = 0)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])