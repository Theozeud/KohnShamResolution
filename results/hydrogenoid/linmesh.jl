using KohnShamResolution

include("../../benchmarktools/Hydrogenoid/setup.jl")


const problinmesh = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 60, 
                            Nmesh         = 70,
                            typemesh      = linmesh, 
                            typebasis     = ShortP1IntLegendreBasis, 
                            optsmesh      = (;),  
                            optsbasis     = (normalize = false, ordermax = 2),                           
                            name          = "IntLeg2-linmesh",
                            nU            = nothing)

const probgeomesh = HydrogenoidProblem(; 
                            T             = Double64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 30, 
                            Nmesh         = 70,
                            typemesh      = geometricmesh, 
                            typebasis     = ShortP1IntLegendreBasis, 
                            optsmesh      = (s = 0.95,),  
                            optsbasis     = (normalize = false, ordermax = 2),                           
                            name          = "IntLeg2-geomesh",
                            nU            = nothing)
#@time λgeo = eigvals_hydro(problinmesh)

#@time sol = eigen_hydro(problinmesh)
##plt = plot_eigenvector(1, Ugeo, probgeomesh)

convNmesh = convergenceNmesh(2 .^(1:8), [problinmesh, probgeomesh]; nums = [1])

convergence_plot_Nmesh(convNmesh)