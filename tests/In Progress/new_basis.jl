using KohnShamResolution
using Test

T = Float64
mesh = linmesh(0,10,100)

@time old_basis = ShortIntLegendreBasis(mesh, T; ordermax = 3)
@time new_basis = IntLegendreGenerator(mesh, T; ordermax = 3)

@time oldM = mass_matrix(old_basis)
@time newM = mass_matrix(new_basis)

@test oldM == newM

@time oldA = stiffness_matrix(old_basis)
@time newA = stiffness_matrix(new_basis)

@test oldA == newA

@time oldM1 = weight_mass_matrix(old_basis,-1)
@time newM1 = weight_mass_matrix(new_basis,-1)

@test oldM1 == newM1

@time oldM2 = weight_mass_matrix(old_basis,-2)
@time newM2 = weight_mass_matrix(new_basis,-2)

@test oldM2 == newM2

@time "oldF" oldF = weight_mass_3tensor(old_basis,Monomial(-1))
@time "newF" newF = weight_mass_3tensor(new_basis,Monomial(-1))

@test oldF == newF
