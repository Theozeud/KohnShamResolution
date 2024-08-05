include("../benchmark tools/include.jl")

# File to save plot
file = "image/hydrogenoide/with short IntLeg - P1 elements/convergence/"

# Basis to test
basis_p1 = ShortP1Basis
basis_il = ShortP1IntLegendreBasis
opts_basis_p1 = (normalize = true, left = false, right = false)                              # Order 1
opts_basis_il2 = (normalize = true, ordermin = 2, ordermax = 2, left = false, right = false) # Order 2
opts_basis_il3 = (normalize = true, ordermin = 2, ordermax = 3, left = false, right = false) # Order 3
opts_basis_il4 = (normalize = true, ordermin = 2, ordermax = 4, left = false, right = false) # Order 4
opts_basis_il5 = (normalize = true, ordermin = 2, ordermax = 5, left = false, right = false) # Order 5
basis = (P1 = basis_p1, IntLeg2 = basis_il, IntLeg3 = basis_il, IntLeg4 = basis_il, IntLeg5 = basis_il)
opts_basis = [opts_basis_p1, opts_basis_il2, opts_basis_il3, opts_basis_il4, opts_basis_il5]

# Formula to get Rmax
# z^1.5 xe^(-zx) = 1e-16
# -zxe^(-zx)= -1e-16/sqrt(z)
# x = -lambertw(-1e-16/sqrt(z), -1)/z

# Parameters to test
# - degree of approximations
# - Cutoff Rmax
# - Parameter l
# - Parameter z

# Plots for l = 0 for different eigenvalues
if true
    T = Float64
    l = 0
    z = 1
    Rmax = Int(round(-lambertw(-1e-16/sqrt(z), -1)/z))
  
    typemesh = linmesh
    vecNmesh     = 2 .^(4:7)

    plt1, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 1)
    #plt2, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 2)
    #plt3, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 3)
    #plt4, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 4)
    #=
    title = "l = $l, z = $z, Rmax = $Rmax, T = $T"
    savefig(plt1, file*"Cvg - "*title*", n = 1")
    savefig(plt2, file*"Cvg - "*title*", n = 2")
    savefig(plt3, file*"Cvg - "*title*", n = 3")
    savefig(plt4, file*"Cvg - "*title*", n = 4")
    =#
end

# Plots for different l for degree 1 to 5 for the the first four eigenvalue
if false
    T = Float64
    z = 1
    Rmax = Int(round(-lambertw(-1e-16/sqrt(z), -1)/z))
    typemesh = linmesh
    vecNmesh     = 2 .^(4:7)
    nb = 1

    plt1, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = 1, nb_eigval = nb)
    plt2, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = 2, nb_eigval = nb)
    plt3, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = 3, nb_eigval = nb)
    plt4, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = 4, nb_eigval = nb)

    title = "n = $nb, z = $z, Rmax = $Rmax, T = $T"
    savefig(plt1, file*"Cvg - "*title*", l = 1")
    savefig(plt2, file*"Cvg - "*title*", l = 2")
    savefig(plt3, file*"Cvg - "*title*", l = 3")
    savefig(plt4, file*"Cvg - "*title*", l = 4")
end

# Plots for different z
if false
    T = Float64
    l = 0
    z = 2
    Rmax = Int(round(-lambertw(-1e-16/sqrt(z), -1)/z))
    typemesh = linmesh
    vecNmesh     = 2 .^(4:7)

    plt1, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 1)
    plt2, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 2)
    plt3, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 3)
    plt4, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 4)

    title = "l = $l, z = $z, Rmax = $Rmax, T = $T"
    savefig(plt1, file*"Cvg - "*title*", n = 1")
    savefig(plt2, file*"Cvg - "*title*", n = 2")
    savefig(plt3, file*"Cvg - "*title*", n = 3")
    savefig(plt4, file*"Cvg - "*title*", n = 4")
end

if false
    T = Float64
    l = 0
    z = 3
    Rmax = Int(round(-lambertw(-1e-16/sqrt(z), -1)/z))
    typemesh = linmesh
    vecNmesh     = 2 .^(4:7)

    plt1, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 1)
    plt2, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 2)
    plt3, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 3)
    plt4, _ = test_convergence_withNmesh(vecNmesh, Rmax, z, basis, typemesh; opts_basis = opts_basis, T = T, l = l, nb_eigval = 4)

    title = "l = $l, z = $z, Rmax = $Rmax, T = $T"
    savefig(plt1, file*"Cvg - "*title*", n = 1")
    savefig(plt2, file*"Cvg - "*title*", n = 2")
    savefig(plt3, file*"Cvg - "*title*", n = 3")
    savefig(plt4, file*"Cvg - "*title*", n = 4")
end