include("../benchmark.jl")

problimesh = hydrogenoid(;  z             = 1, 
                            l             = 0, 
                            Rmax          = 40, 
                            Nmesh         = 60,
                            typemesh      = linmesh, 
                            typebasis     = ShortP1IntLegendreBasis, 
                            optsmesh      = NamedTuple(), 
                            optsbasis     = (normalize = true, ordermax = 4, ordermin = 2), 
                            T             = Float64, 
                            name          = "IntLeg4")

Ulin, λlin = eigen_hydro(problimesh)

probgeomesh = hydrogenoid(; z             = 1, 
                            l             = 0, 
                            Rmax          = 40, 
                            Nmesh         = 60,
                            typemesh      = geometricmesh, 
                            typebasis     = ShortP1IntLegendreBasis, 
                            optsmesh      = (s = 0.9,), 
                            optsbasis     = (normalize = false, ordermax = 4, ordermin = 2), 
                            T             = Float64, 
                            name          = "IntLeg2")

Ugeo, λgeo = eigen_hydro(probgeomesh)

vecNmesh = 2 .^ [2,3,4,5,6,7,8]
println("Calcul erreur de convergence pour grille uniforme")
Errlin = convergenceNmesh(vecNmesh, problimesh; num_eig = 1)

println("Calcul erreur de convergence pour grille géométrique")
Errgeo = convergenceNmesh(vecNmesh, probgeomesh; num_eig = 1)

plt = plot( size = (800,600), margin = 0.5Plots.cm, legend = :bottomleft, xaxis=:log, yaxis=:log,
            legendfontsize  = 12,  
            titlefontsize   = 12,
            guidefontsize   = 12,
            tickfontsize    = 12)
xlabel!(plt, "Nmesh")
ylabel!(plt, "Error on the 1-th eigenvalues")
title!(plt, "Rmax = 40, Nmesh = 60, IntLeg4, z = 1, l = 0")
plot!(plt, vecNmesh, Errlin, lw = 4, label = "linmesh", markershape = :x, markersize = 10)
plot!(plt, vecNmesh, Errgeo, lw = 4, label = "geomesh", markershape = :circle, markersize = 10)
plt