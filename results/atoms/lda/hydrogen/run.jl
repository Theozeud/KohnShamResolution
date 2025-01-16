include("../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig()

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.7), 
                model           = LSDA(1, 1), 
                Rmax            = 40.0, 
                Nmesh           = 40,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 3,),
                typediscre      = LSDADiscretization, 
                name            = "test", 
                scftol          = 1e-6,
                maxiter         = 1,
                hartree         = true,
                logconfig       = logconfig)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])