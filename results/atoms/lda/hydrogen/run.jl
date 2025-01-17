include("../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.5), 
                model           = LSDA(1, 1), 
                Rmax            = 30.0, 
                Nmesh           = 50,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 3,),
                typediscre      = LSDADiscretization, 
                name            = "test", 
                scftol          = 1e-6,
                maxiter         = 50,
                hartree         = true,
                logconfig       = logconfig)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])