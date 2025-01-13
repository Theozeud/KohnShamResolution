include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig()

problem = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.7), 
                model           = ReducedHartreeFock(1, 1), 
                Rmax            = 40.0, 
                Nmesh           = 60,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 3,), 
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])