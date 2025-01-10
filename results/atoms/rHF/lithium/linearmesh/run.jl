include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# MODEL
z = 3
N = 3

# LOG CONFIG
logconfig = LogConfig()

problem = AtomProblem(;
                T               = Float64, 
                lh              = 4, 
                method          = CDA(0.7), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 60.0, 
                Nmesh           = 60,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = ShortP1IntLegendreBasis, 
                optsbasis       = (ordermax = 4,), 
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])