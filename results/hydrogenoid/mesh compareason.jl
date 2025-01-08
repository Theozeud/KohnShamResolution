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
                            optsbasis     = (normalize = false, ordermax = 4),                           
                            name          = "IntLeg2-linmesh",
                            nU            = nothing)

const probgeomesh = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 60, 
                            Nmesh         = 70,
                            typemesh      = geometricmesh, 
                            typebasis     = ShortP1IntLegendreBasis, 
                            optsmesh      = (s = 0.9,),  
                            optsbasis     = (normalize = false, ordermax = 4),                           
                            name          = "IntLeg2-geomesh",
                            nU            = nothing)
#@time λgeo = eigvals_hydro(problinmesh)

#@time sol = eigen_hydro(problinmesh)
##plt = plot_eigenvector(1, Ugeo, probgeomesh)

convNmesh = convergenceNmesh(2 .^(3:8), [problinmesh, probgeomesh]; nums = [4])

convergence_plot_Nmesh(convNmesh)