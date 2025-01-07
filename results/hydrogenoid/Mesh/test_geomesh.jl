include("../benchmark.jl")

probgeomesh = hydrogenoid(; z             = 1, 
                            l             = 0, 
                            Rmax          = 60, 
                            Nmesh         = 70,
                            typemesh      = geometricmesh, 
                            typebasis     = ShortP1IntLegendreBasis, 
                            optsmesh      = (s = 0.9,), 
                            optsbasis     = (normalize = false, ordermax = 4, ordermin = 2), 
                            T             = Double64, 
                            name          = "IntLeg4")

@time λgeo = eigvals_hydro(probgeomesh)

#Ugeo, λgeo = eigen_hydro(probgeomesh)
#plt = plot_eigenvector(1, Ugeo, probgeomesh)

λgeo