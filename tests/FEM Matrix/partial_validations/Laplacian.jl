using KohnShamResolution
using LinearAlgebra
using Plots
using UnPack
using GenericLinearAlgebra
using DoubleFloats

# Structure Problem
struct LaplacianProblem
    name
    T
    typebasis
    optsbasis
    typemesh
    optsmesh
    Rmax
    Nmesh
end

# Structure Solution
struct LaplacianSolution
    prob
    λ
    u
end

# Theoretical Eigenvalue and Eigenvectors
theoretical_eigenvector(n, Rmax, T) = x-> sqrt(T(2)/T(Rmax)) * sin(n*T(π)*x/T(Rmax)) 
theoretical_eigenvalue(n, Rmax, T)  = (T(π)/T(Rmax))^2 * n^2

# Compute eigenvectors and eigenvalues
function eigen_(problems)
    sols = LaplacianSolution[]
    for problem ∈ problems
        @unpack T, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis = problem
        mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
        basis = typebasis(mesh, T; optsbasis...)
        M₀  = Symmetric(mass_matrix(basis))
        A   = Symmetric(mass_matrix(deriv(basis)))
        u, λ = eigen(A, M₀)
        push!(sols, LaplacianSolution(problem, λ, u))
    end
    sols
end

# Compute eigen values (usefull for Double64 for which eigen is not defined)
function eigvals_(problems)
    sols = LaplacianSolution[]
    for problem ∈ problems
        @unpack T, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis = problem
        mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
        basis = typebasis(mesh, T; optsbasis...)
        M₀  = mass_matrix(basis)
        A   = Symmetric(mass_matrix(deriv(basis)))
        λ = eigvals(inv(M₀) * A)
        push!(sols, LaplacianSolution(problem, λ, nothing))
    end
    sols
end

# Plot Eigenvalue
function plot_eigenvalue(sols)
    plt = plot( size= (1500, 700), margin = 0.5Plots.cm,                
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
    mindex = 0
    RMAX   = 0
    for sol ∈ sols
        plot!(plt, sol.λ, label = sol.prob.name, lw = 4)
        if length(sol.λ) > mindex
            mindex = length(sol.λ)
        end
        if sol.prob.Rmax > RMAX
            RMAX = sol.prob.Rmax
        end
    end
    plot!(plt, theoretical_eigenvalue.(1:mindex, RMAX, sol.prob.T), label = "Theoretical", lw = 4)
    return plt
end

##############################################################################
# Plot Eigenvector

function plot_eigenvector(sols, n; resol = 5000)
    plt = plot( size = (1500, 700), margin = 1Plots.cm,                
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
    xlabel!(plt, "r")
    title!(plt, string(n)*"-th eigenvector")
    RMAX = 0
    for sol ∈ sols
        eign = build_on_basis(basis, sol.u[:,n])
        eign_n = eign / sqrt(integrate(eign*eign, zero(T), T(Rmax)))
        X = LinRange(O, sol.prob.Rmax, resol)
        plot!(plt, X, (sign(theoretical_eigenvector(n, Rmax, sol.prob.T)(1) * eign_n(1)) * eign_n).(X), label = sol.prob.name, lw = 4)
        if sol.prob.Rmax > RMAX
            RMAX = sol.prob.Rmax
        end
    end
    X = LinRange(O, RMAX, resol)
    plot!(plt, X, theoretical_eigenvector(n, Rmax, sol.prob.T).(X), label = "Theoretical n = $n", lw = 4)
    plt
end

##################################################################################
# Compute Error

function compute_error_eigenvalue(vecNmesh, Rmax, num, T)

    label = ["Différence finie", "P1", "IntLeg 2", "IntLeg 3", "IntLeg 4"]
    ϵerror = zeros(T, length(vecNmesh), length(label))
    for (i, Nmesh) ∈ enumerate(vecNmesh)
        # Creation of the mesh
        m = linmesh(Rmin, Rmax, Nmesh, T)
        
        # True eigenvalue
        @time "Nmesh = "*string(Nmesh) begin
        true_val = val_theo.(1:Nmesh, Rmax, T)

        # With finite difference
        vals_diff = eigen_with_finite_diff(Nmesh, Rmax)

        # With P1
        vals_p1 = eigen_with_P1(m, T)

        # With P1-Integrated Legendre Polynomials ordre 2
        vals_il2 = eigen_with_P1IntLeg(m, 2, T)

        # With P1-Integrated Legendre Polynomials ordre 3
        vals_il3 = eigen_with_P1IntLeg(m, 3, T)

        # With P1-Integrated Legendre Polynomials ordre 4
        vals_il4 = eigen_with_P1IntLeg(m, 4, T)
        end

        # Compute the error for eigenvalues and the fundamental
        ϵerror[i,1] = abs(vals_diff[num] - true_val[num])
        ϵerror[i,2] = abs(vals_p1[num] - true_val[num])
        ϵerror[i,3] = abs(vals_il2[num] - true_val[num])
        ϵerror[i,4] = abs(vals_il3[num] - true_val[num])
        ϵerror[i,5] = abs(vals_il4[num] - true_val[num])
    end

    # Creation of the plot for eigenvalue
    plt_ϵ_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
                        legendfontsize  = 14,  
                        titlefontsize   = 14,
                        guidefontsize   = 14,
                        tickfontsize    = 14)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the $num-th eigenvalue")
    title!(plt_ϵ_error, "Rmax = $Rmax")
    for i ∈ eachindex(label)
        plot!(plt_ϵ_error, vecNmesh, ϵerror[:, i], lw = 4, label = label[i], markershape = :x, markersize = 10)
    end
    plt_ϵ_error
end

##############################################
prob = LaplacianProblem("problem", 
                        Float64,
                        ShortP1IntLegendreBasis,
                        (normalize = false, ordermax = 4),
                        linmesh,
                        (;),
                        10,
                        100)

sol = eigvals_([prob])

plot_eigenvalue(sol)

#=
plt_eigenvalue = plot_eigenvalue(100, 10, [2,3])

savefig(plt_eigenvalue, "test/image_tests/Laplacien - Valeurs Propres 1-2-3")
plt_eigenvalue = plot_eigenvalue(100,10, [2,3,4,5])
savefig(plt_eigenvalue, "test/image_tests/Laplacien - Valeurs Propres 1-2-3-4-5")


plt_eigenvector = plot_eigenvector(100,10, 198, 3)
savefig(plt_eigenvector, "test/image_tests/Laplacien - Vecteur Propre 198 - ordre 3")
plt_eigenvector = plot_eigenvector(100,10, 199, 3)
savefig(plt_eigenvector, "test/image_tests/Laplacien - Vecteur Propre 199 - ordre 3")
plt_eigenvector = plot_eigenvector(100,10, 200, 3)
savefig(plt_eigenvector, "test/image_tests/Laplacien - Vecteur Propre 200 - ordre 3")

plt_eigenvalue_error = compute_error_eigenvalue(vecNmesh, 10, 1)
savefig(plt_eigenvalue_error, "test/image_tests/Laplacien - Erreurs valeurs propres")

compute_error_eigenvalue(vecNmesh, 10, 1)
=#
