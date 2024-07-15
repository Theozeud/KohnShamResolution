using KohnShamResolution
using LinearAlgebra
using Plots
using GenericLinearAlgebra
using DoubleFloats

# Theoretical Eigenvalue and Eigenvectors
function eigval_theo(T, n, z)
    -T(z)^2/(T(2)*T(n)^2)
end

##############################################################################
# Function to compute eigenvalue and eigenvector for different discretizations

# With P1
function eigen_with_P1(T, mesh)
    normalize = true
    left = false
    right = false
    basis = ShortP1Basis(mesh, T; left = left, right = right, normalize = normalize)
    deriv_basis = deriv(basis)
    A   = Symmetric(mass_matrix(deriv_basis))
    M₀  = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    H = T(0.5) * A  -  M₋₁
    eigvals(inv(M₀) * H)
end

# With P1-Integrated Legendre Polynomials
function eigen_with_P1IntLeg(T, mesh, ordermax)
    normalize = true
    ordermin = 2
    left = false
    right = false
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize, left = left, right = right)
    deriv_basis = deriv(basis)
    A   = Symmetric(mass_matrix(deriv_basis))
    M₀  = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    H = T(0.5) * A  -  M₋₁
    eigvals(inv(M₀) * H)
end

# With given basis
function eigvals_hydro(basis)
    deriv_basis = deriv(basis)
    A   = Symmetric(mass_matrix(deriv_basis))
    M₀  = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    H = T(0.5) * A  -  M₋₁
    eigvals(inv(M₀) * H)
end

##################################################################################
# Compute Error

function test_convergence_withNmesh(vecNmesh::AbstractVector, Rmax::Real, vecBasis::NamedTuple, typemesh; opts_mesh = NamedTuple(), opts_basis = [NamedTuple() for i ∈ eachindex(vecBasis)], T = Float64, nb_eigval = 1)
    plt_ϵ_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :outertopright, xaxis=:log, yaxis=:log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the "*string(nb_eigval)*"-th eigenvalues")
    title!(plt_ϵ_error, "Rmax = "*string(Rmax)*" and z = 1")
    ϵerror = zeros(T, length(vecBasis), length(vecNmesh))
    for (i, Basis) ∈ enumerate(vecBasis)
        println("Basis "*string(keys(vecBasis)[i]))
        (_, ϵ_error) = test_convergence_withNmesh(vecNmesh, Rmax, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis[i], T = T, nb_eigval = nb_eigval)
        plot!(plt_ϵ_error, vecNmesh, ϵ_error, lw = 3, label = "Basis "*string(keys(vecBasis)[i])*", ϵ"*string(nb_eigval), markershape = :x, markersize = 10)
        ϵerror[i, : ] = ϵ_error
    end
    (plt_ϵ_error, ϵerror)
end

function test_convergence_withNmesh(vecNmesh::AbstractVector, Rmax::Real, Basis, vectypemesh::NamedTuple; opts_mesh = [NamedTuple() for i ∈ eachindex(vectypemesh)], opts_basis = NamedTuple(), T = Float64, nb_eigval = 1)
    plt_ϵ_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :outertopright, xaxis=:log, yaxis=:log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the "*string(nb_eigval)*"-th eigenvalues")
    title!(plt_ϵ_error, "Rmax = "*string(Rmax)*" and z = 1")
    ϵerror = zeros(T, length(vectypemesh), length(vecNmesh))
    for (i, typemesh) ∈ enumerate(vectypemesh)
        println("mesh "*string(keys(vectypemesh)[i]))
        (_, ϵ_error) = test_convergence_withNmesh(vecNmesh, Rmax, Basis, typemesh; opts_mesh = opts_mesh[i], opts_basis = opts_basis, T = T, nb_eigval = nb_eigval)
        plot!(plt_ϵ_error, vecNmesh, ϵ_error, lw = 3, label = "mesh "*string(keys(vectypemesh)[i])*", ϵ"*string(nb_eigval), markershape = :x, markersize = 10)
        ϵerror[i, : ] = ϵ_error
    end
    (plt_ϵ_error, ϵerror)
end


function test_convergence_withNmesh(vecNmesh, Rmax::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, nb_eigval = 1)
    Rmin = zero(T)
    ϵ_error = zeros(T, length(vecNmesh))
    for (i,Nmesh) ∈ enumerate(vecNmesh)
        m = typemesh(Rmin, Rmax, Nmesh; T = T, opts_mesh...)
        basis = Basis(m, T; opts_basis...)
        @time "Nmesh = $Nmesh" ϵ = eigvals_hydro(basis)
        ϵ_error[i] = abs(ϵ[nb_eigval] - eigval_theo(T, nb_eigval, 1))
    end
    # Creation of the plot for eigenvalue
    plt_ϵ_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :outertopright, xaxis=:log, yaxis=:log,
                                 legendfontsize  = 12,  
                                 titlefontsize   = 12,
                                 guidefontsize   = 12,
                                 tickfontsize    = 12)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the $nb_eigval-th eigenvalues")
    title!(plt_ϵ_error, "Rmax = $Rmax and z = 1")
    plot!(plt_ϵ_error, vecNmesh, ϵ_error, lw = 3, label = "ϵ$nb_eigval", markershape = :x, markersize = 10)
    return (plt_ϵ_error, ϵ_error)
end




#=
# Creation of the model
z = 1

# Discretization 
Rmin = 0
Rmax = 100
Nmesh = 100
m = logmesh(Rmin, Rmax, Nmesh)
basis = ShortP1IntLegendreBasis(m; left = false, right = false, normalize = true, ordermin = 2, ordermax = 3)

# Solve
deriv_basis = deriv(basis)
A   = mass_matrix(deriv_basis)
M₀  = mass_matrix(basis)
M₋₁ = weight_mass_matrix(basis, -1)


H = Symmetric(1/2 * A  - z .* M₋₁)
@time ϵ, U = eigen(H,M₀)
=#
