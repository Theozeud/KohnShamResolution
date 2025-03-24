using KohnShamResolution

include("../../benchmarktools/hydrogenoid/setup.jl")


const problinmesh = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 90, 
                            Nmesh         = 70,
                            typemesh      = linmesh, 
                            typebasis     = P1IntLegendreGenerator, 
                            optsmesh      = (),  
                            optsbasis     = (ordermax = 2, ),                           
                            name          = "IntLeg2-linmesh",
                            nU            = nothing)

const probgeomesh = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 90, 
                            Nmesh         = 70,
                            typemesh      = geometricmesh, 
                            typebasis     = P1IntLegendreGenerator, 
                            optsmesh      = (s = 0.9,),  
                            optsbasis     = (ordermax = 2, ),                           
                            name          = "IntLeg2-geomesh",
                            nU            = nothing)

const probexpmesh = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 90, 
                            Nmesh         = 70,
                            typemesh      = expmesh, 
                            typebasis     = P1IntLegendreGenerator, 
                            optsmesh      = (s = 2.0,),  
                            optsbasis     = (ordermax = 2, ),                           
                            name          = "IntLeg2-expmesh",
                            nU            = nothing)
#@time λgeo = eigvals_hydro(problinmesh)

#@time sol = eigen_hydro(problinmesh)
##plt = plot_eigenvector(1, Ugeo, probgeomesh)

convNmesh = convergenceNmesh(2 .^(3:10), [problinmesh, probgeomesh, probexpmesh]; nums = [1])

convergence_plot_Nmesh(convNmesh)