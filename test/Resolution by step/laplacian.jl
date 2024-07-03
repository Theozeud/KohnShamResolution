using KohnShamResolution
using Plots
using LinearAlgebra

# General Discretization Parameters
T = Float64
Rmin = 0
Rmax = 1000
Nmesh = 500
m = linmesh(Rmin,Rmax,Nmesh)
normalize = true

# Theoretical Eigenvalue and Eigenvectors
vect_theo(n) = x-> sqrt(2/Rmax) * sin(n*π*x/Rmax) 
val_theo(n)  = (π/Rmax)^2 * n^2


# With finite difference
Nmesh = 100
A = SymTridiagonal(2*ones(Nmesh), -ones(Nmesh-1)) .* ((Nmesh - 1)/ Rmax)^2
vals_diff, eigs_diff = eigen(A)


# With P1
left = false
right = false
basis = ShortP1Basis(m, T; left = left, right = right, normalize = normalize)
AP1   = mass_matrix(deriv(basis))
vals_p1, eigs_p1 = eigen(AP1)


# With P1-Integrated Legendre Polynomials
ordermin = 2
ordermax = 3
basis = ShortP1IntLegendreBasis(m, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize)
AIntleg   = mass_matrix(deriv(basis))
vals_il, eigs_il = eigen(AIntleg)

# Plot Eigenvalue
plt_eigenvalue = plot(size= (1000, 500))
plot!(vals_diff[1:Nmesh], label = "Finite Difference", lw = 2)
plot!(vals_p1[1:Nmesh], label = "P1", lw = 2)
plot!(vals_il[1:Nmesh], label = "P1 + Integrated Legendre", lw = 2)
plot!(val_theo.(1:Nmesh), label = "Theoretical", lw = 2)

plt_eigenvalue_intleg = plot(size= (1000, 500))
plot!(vals_il, label = "P1 + Integrated Legendre", lw = 2)
plot!(val_theo, label = "Theoretical", lw = 2)

plot(plt_eigenvalue, plt_eigenvalue_intleg, layout = (1,2))


#= To build eigenvector on finite element basis
X = LinRange(Rmin, Rmax, Nmesh * 100)
plt_vectA = plot()
for i ∈ 1:min(5, length(basis))
    eigi = build_on_basis(basis, eigs[:,i])
    eigi = eigi / sqrt(integrate(eigi*eigi, 0, Rmax))
    @show integrate(eigi*eigi, 0, Rmax)
    plot!(plt_vectA, X, eigi.(X), label = "num i = "*string(i))
    plot!(plt_vectA, X, vect_theo(i).(X), label = "theo i = "*string(i))
end
plt_vectA
=#


