include("../benchmark tools/include.jl")

# File to save plot

# Fixed Parameters
lₕ = 0
T = Float64
typemesh = logmesh
basis = ShortP1IntLegendreBasis
ordermax = 3
opts_basis = (normalize = true, ordermin = 2, ordermax = ordermax, left = false, right = false)

# Default parameters
Rmax = 100
Nmesh = 40
z = 1

# Study with respect to
vecz     = 1:1:4

# Plots
plt, _ = hydrogenoide_test_eigenvalue(Rmax, Nmesh, vecz, basis, typemesh; opts_basis = opts_basis, T = T, lₕ = lₕ)

savefig(plt, "image/hydrogenoide/with short IntLeg - P1 elements/Valeurs propres avec ordremax = "*string(ordermax)*", Nmesh = "*string(Nmesh)*" et Rmax = "*string(Rmax))