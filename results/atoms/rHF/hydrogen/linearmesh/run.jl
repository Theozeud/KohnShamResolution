include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.7), 
                model           = ReducedHartreeFock(1, 1), 
                Rmax            = 80, 
                Nmesh           = 50,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 8,
                hartree         = true,
                degen_tol       = 1e-15,
                logconfig       = logconfig)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])