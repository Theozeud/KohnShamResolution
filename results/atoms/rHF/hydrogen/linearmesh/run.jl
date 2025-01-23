include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig()

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.5), 
                model           = SlaterXα(1, 1), 
                Rmax            = 50, 
                Nmesh           = 70,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 3,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-12,
                maxiter         = 150,
                hartree         = true,
                logconfig       = logconfig)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])