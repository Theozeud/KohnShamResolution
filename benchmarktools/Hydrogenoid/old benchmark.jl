using KohnShamResolution
using LinearAlgebra
using Plots
using GenericLinearAlgebra
using DoubleFloats
using UnPack

#############################################################################################
#                                    STRUCTURE HYDROGENOID
#############################################################################################
struct Hydrogenoid
    z
    l
    Rmax
    Nmesh
    mesh
    basis
    name
    conv 
end

function hydrogenoid(; z, l, Rmax, Nmesh, typemesh, typebasis, optsmesh, optsbasis, T, name)
    mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
    basis = typebasis(mesh, T; optsbasis...)
    conv(N) = hydrogenoid(;z = z, l = l, Rmax = Rmax, Nmesh = N, typemesh = typemesh, typebasis = typebasis, optsmesh = optsmesh, optsbasis = optsbasis, T = T, name = name)
    Hydrogenoid(z, l, Rmax, Nmesh, mesh, basis, name, conv)
end

function hydrogenoid_conv(; z, l, Rmax, typemesh, typebasis, optsmesh, optsbasis, T, name)
    N -> hydrogenoid(;z = z, l = l, Rmax = Rmax, Nmesh = N, typemesh = typemesh, typebasis = typebasis, optsmesh = optsmesh, optsbasis = optsbasis, T = T, name = name)
end

#############################################################################################
#                                     THEORETICAL DATA
#############################################################################################
# Theoretical Eigenvalues
function eigval_theo(n, z, l, T = Float64)
    -T(z)^2/(T(2)*T(n+l)^2)
end

# Theoretical Eigenvectors
function eigvect_theo(n, z, l, T = Float64)
    exp_part(x) = exp(-z*abs(x))*x
    C = zeros(T, n-l)
    C[1] = z^(3/2)/sqrt(π)
    for i ∈ 2:n-l
        C[i]= -2/n * (n-l-i)/(i*(2*l+i+1)) * C[i-1]
    end
    return x -> begin
        val = C[n-l] 
        for j ∈ n-l-1:-1:1
            val = val * x + C[j]
        end
        val * exp_part(x) * x^l
    end
end

#############################################################################################
#                                          RESOLUTION
#############################################################################################

# Normamization procedure
function L2normalization(U, M₀)
    ax = axes(U, 1)
    nU = zero(U)
    for k ∈ axes(U,2)
        normalization = sqrt(sum([U[i,k] * U[j,k] * M₀[i,j] for i∈ax for j∈ax]) * 4π)
        nU[:,k] .= U[:,k] / normalization
    end
    nU
end

# Eigenvalues and Eigenvectors with IntLeg
function eigen_hydro(hydro)
    @unpack z, l, mesh, basis = hydro
    T = KohnShamResolution.bottom_type(basis)
    deriv_basis = deriv(basis)
    A   = Symmetric(mass_matrix(deriv_basis))
    M₀  = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    if l == 0
        H = T(0.5) * A  -  z*M₋₁
        Λ, U = eigen(H, M₀)
        nU = L2normalization(U, M₀)
        return  (nU, Λ)
    else
        M₋₂ = Symmetric(weight_mass_matrix(basis, -2))
        H = T(0.5) * A  - z * M₋₁ + T(0.5) * l*(l+1) * M₋₂
        Λ, U = eigen(H, M₀)
        nU = L2normalization(U, M₀)
        return  (nU, Λ)
    end
end

# Eigenvalue with Eigvals (usefull for Double64)
function eigvals_hydro(hydro)
    @unpack z, l, mesh, basis = hydro
    T = KohnShamResolution.bottom_type(basis)
    deriv_basis = deriv(basis)
    A   = Symmetric(mass_matrix(deriv_basis))
    M₀  = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    if l == 0
        H = T(0.5) * A  -  z * M₋₁
        return  eigvals(inv(M₀) * H)
    else
        M₋₂ = Symmetric(weight_mass_matrix(basis, -2))
        H = T(0.5) * A  -  z * M₋₁ + T(0.5) * l*(l+1) * M₋₂
        return  eigvals(inv(M₀) * H)
    end
end


#############################################################################################
#                                     PLOT EIGENVECTOR
#############################################################################################
function eval_on_basis(U, basis, x)
    T = bottom_type(basis)
    ux = zero(T)
    for i ∈ eachindex(U)
        ux += U[i] * KohnShamResolution.eval_basis(basis, i, x)
    end
    ux
end

function plot_eigenvector(n, U, hydro)
    @unpack z, l, mesh, basis, name = hydro
    Rmin = first(mesh.points)
    Rmax = last(mesh.points)
    Nmesh = length(mesh.points)
    T = bottom_type(basis)
    plt = plot( size = (800,600), margin = 0.5Plots.cm, legend = :topright,
                legendfontsize  = 13,  
                titlefontsize   = 13,
                guidefontsize   = 13,
                tickfontsize    = 13)
    xlabel!(plt, "r")
    ylabel!(plt, "$n-th eigenvector")
    title!(plt, "Rmax = $Rmax, Nmesh = $Nmesh, $name, z = $z, l = $l")
    X = LinRange(Rmin, Rmax, 1000)
    
    solnum = [sign(eigvect_theo(n, z, l, T)(x)) *abs(eval_on_basis(U[:,n], basis, x)) for x ∈ X]
    plot!(X, solnum, label = "Numérique", lw = 3)
    plot!(X, eigvect_theo(n, z, l, T).(X), lw = 2, label = "Théorique", ls = :dash, color = :black)
end

#############################################################################################
#                                     CONVERGENCE WITH NMESH
#############################################################################################

function convergenceNmesh(vecNmesh::AbstractVector, hydro; num_eig = 1)
    @unpack basis, z, l, Rmax, name = hydro
    T = bottom_type(basis)
    Err = zeros(T, length(vecNmesh))
    for (i, Nmesh) ∈ enumerate(vecNmesh)
        @time "Nmesh = $Nmesh" U, λ = eigen_hydro(hydro.conv(Nmesh))
        Err[i] = abs(λ[num_eig] - eigval_theo(num_eig, z, l, T))
    end
    title = "Rmax = $Rmax, $name, z = $z, l = $l"
    Err, title
end