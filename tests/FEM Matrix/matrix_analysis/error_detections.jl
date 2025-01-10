using KohnShamResolution
using LinearAlgebra

T = Float64

Rmin = 0
Rmax  = 30
Nmesh = 10
mesh = geometricmesh(Rmin, Rmax, Nmesh; T = T, s = 0.7)


function meshsize(mesh::Mesh)
    mesh.points[2:end] - mesh.points[1:end-1]
end

basis = ShortIntLegendreBasis(mesh, T; ordermax = 2)

M₀  = mass_matrix(basis)
@show any(isnan, M₀)
@show any(isinf, M₀)
A   = stiffness_matrix(basis)
@show any(isnan, A)
@show any(isinf, A)
M₋₁ = weight_mass_matrix(basis, -1)
@show any(isnan, M₋₁)
@show any(isinf, M₋₁)
M₋₂ = weight_mass_matrix(basis, -2)
@show any(isnan, M₋₂)
@show any(isinf, M₋₂)

