using KohnShamResolution
using LinearAlgebra
using GenericLinearAlgebra
using DoubleFloats

N = 200

T = Float64
Rmin = T(0)
Rmax = T(50)
m = geometricmesh(Rmin, Rmax, N; T = T, s = 0.9)
Fbasis = ShortP1IntLegendreBasis(m, T; normalize = false, ordermax = 5)
deriv_basis = deriv(Fbasis)
@time "FA" FA   = mass_matrix(deriv_basis)
@time "FM₀" FM₀  = mass_matrix(Fbasis)
@time "FM₋₁" FM₋₁ = weight_mass_matrix(Fbasis, -1)
@time "FM₋₂" FM₋₂ = weight_mass_matrix(Fbasis, -2)
@time "FF" FF = weight_mass_3tensor(Fbasis, Monomial(-1))

T = Double64
Rmin = T(0)
Rmax = T(50)
m = geometricmesh(Rmin, Rmax, N; T = T, s = 0.9)
Dbasis = ShortP1IntLegendreBasis(m, T; normalize = false, ordermax = 5)
deriv_basis = deriv(Dbasis)
@time "DA" DA   = mass_matrix(deriv_basis)
@time "DM₀" DM₀  = mass_matrix(Dbasis)
@time "DM₋₁" DM₋₁ = weight_mass_matrix(Dbasis, -1)
@time "DM₋₂" DM₋₂ = weight_mass_matrix(Dbasis, -2)
@time "DF" DF = weight_mass_3tensor(Dbasis, Monomial(-1))


@show norm(DA - FA)
@show norm(DM₀ - FM₀)
@show norm(DM₋₁ - FM₋₁)
@show norm(DM₋₂ - FM₋₂)
@show norm(DF  - FF)