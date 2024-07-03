include("../benchmark tools/include.jl")

# File to save plot

# Fixed Parameters
T = Float64
typemesh = logmesh
basis = ShortIntLegendreBasis
ordermax = 3
opts_basis = (normalize = true, ordermin = 1, ordermax = ordermax)

# Default parameters
Rmax = 100
Nmesh = 40
z = 1

# Study with respect to
vecRmax  = 20:5:150
vecNmesh = 10:10:200
vecz     = 1:1:20

# Plots
plt_Rmax, _  = conditionning(vecRmax, Nmesh, z, basis, typemesh; opts_basis = opts_basis, T = T)
plt_Nmesh, _ = conditionning(Rmax, vecNmesh, z, basis, typemesh; opts_basis = opts_basis, T = T)
plt_z, _     = conditionning(Rmax, Nmesh, vecz, basis, typemesh; opts_basis = opts_basis, T = T)

## Save Plots

pltfin = plot(plt_Rmax, plt_Nmesh, plt_z, layout = (2,2), size = (2400,1000), margin = 1.2Plots.cm)

savefig(pltfin, "image/hydrogenoide/with short IntLeg elements/Conditionnement pour ordremax = "*string(ordermax))