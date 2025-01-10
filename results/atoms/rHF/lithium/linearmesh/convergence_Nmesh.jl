include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# MODEL
z = 3
N = 3
model = ReducedHartreeFock(z, N)

# LOG CONFIG
logconfig = LogConfig(energy = false)

# CONVERGENCE WITH RESPECT TO NMESH

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.7), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 60.0, 
                Nmesh           = 50,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = ShortP1IntLegendreBasis, 
                optsbasis       = (ordermax = 4,), 
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)
                        
error = convergenceNmesh(2 .^(2:7), [problem])

convergencePlotNmesh(error)



