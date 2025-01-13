include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# MODEL
z = 8
N = 8
model = ReducedHartreeFock(z, N)

# LOG CONFIG
logconfig = LogConfig(energy = false)

# CONVERGENCE WITH RESPECT TO NMESH

problem1 = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.2), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 100.0, 
                Nmesh           = 50,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "0.2", 
                scftol          = 1e-9,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

problem2 = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.5), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 100.0, 
                Nmesh           = 50,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "0.5", 
                scftol          = 1e-9,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

problem3 = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.7), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 100.0, 
                Nmesh           = 50,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "0.7", 
                scftol          = 1e-9,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)
                        
error = convergenceNmesh(2 .^(2:7), [problem1,problem2,problem3]; nums = [1])

convergencePlotNmesh(error)



