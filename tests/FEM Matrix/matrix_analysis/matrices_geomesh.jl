using KohnShamResolution
using LinearAlgebra

T = Float64

Rmin = 0
Rmax  = 30
Nmesh = 200
mesh = geometricmesh(Rmin, Rmax, Nmesh; T = T, s = 0.9)

basis = ShortP1IntLegendreBasis(mesh, T; ordermax = 4)

A = stiffness_matrix(basis)