include("../benchmark tools/include.jl")

# File to save plot

# Fixed Parameters
lₕ = 0
T = Float64
typemesh = logmesh
basis = ShortP1IntLegendreBasis
ordermax = 2
opts_basis = (normalize = true, ordermin = 2, ordermax = ordermax, left = false, right = false)
#opts_mesh = (z = 1)

basis2 = ShortP1Basis
opts_basis2 = (normalize = true, left = false, right = false)

# Default parameters
Rmax = 10
z = 1

# Study with respect to
vecNmesh     = 2 .^(2:9)

# Plots
#plt_ϵ, plt_u, ϵerror, uerror = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, lₕ = lₕ)
plt_ϵ, plt_u, ϵerror, uerror = test_convergence_withNmesh(vecNmesh, Rmax, z, (IntLeg = basis, P1 = basis2), typemesh; opts_basis = [opts_basis, opts_basis2], T = T, lₕ = lₕ)

#savefig(plt, "image/hydrogenoide/with short IntLeg - P1 elements/Convergence pour ordremax = "*string(ordermax)*", z = "*string(z)*" et Rmax = "*string(Rmax))

plt_ϵ