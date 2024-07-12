include("../benchmark tools/include.jl")

# File to save plot

# Fixed Parameters
T = Float64
typemesh = logmesh
basis = ShortP1IntLegendreBasis
ordermax = 4
opts_basis = (normalize = true, ordermin = 2, ordermax = ordermax, left = false, right = false)

# Default parameters
Rmax = 100
Nmesh = 256
z = 1

# Study with respect to
vecRmax  = 20:5:300
vecNmesh = 10:5:400
vecz     = 1:1:20

# Plots
plt_Rmax, _  = conditionning(vecRmax, Nmesh, z, basis, typemesh; opts_basis = opts_basis, T = T)
plt_Nmesh, condNmesh = conditionning(Rmax, vecNmesh, z, basis, typemesh; opts_basis = opts_basis, T = T)
plt_z, _     = conditionning(Rmax, Nmesh, vecz, basis, typemesh; opts_basis = opts_basis, T = T)

## Save Plots

pltfin = plot(plt_Rmax, plt_Nmesh, plt_z, layout = (3,1), size = (1200, 1500), margin = 1.2Plots.cm)

savefig(pltfin, "image/hydrogenoide/with short IntLeg - P1 elements/Conditionnement pour ordremax = $ordermax")