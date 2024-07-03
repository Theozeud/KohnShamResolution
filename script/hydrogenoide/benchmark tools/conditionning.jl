# This files provides function to compute the conditionning of each matrix involved in solving the Hydrogenoid Schrodinger equation

function conditionning(vecRmax::AbstractVector, Nmesh::Int, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64)
    Cond = zeros(T, length(vecRmax), 4)
    for (i,Rmax) ∈ enumerate(vecRmax)
        Cond[i, :] = conditionning(Rmax, Nmesh, z, Basis, typemesh; opts_mesh, opts_basis, T)
    end
    plt = plot_conditionning(Cond, vecRmax, "Rmax", "z = "*string(z)*", Nmesh = "*string(Nmesh))
    (plt, Cond)
end

function conditionning(Rmax::Real, vecNmesh::AbstractVector, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64)
    Cond = zeros(T, length(vecNmesh), 4)
    for (i, Nmesh) ∈ enumerate(vecNmesh)
        Cond[i, :] = conditionning(Rmax, Nmesh, z, Basis, typemesh; opts_mesh, opts_basis, T)
    end
    plt = plot_conditionning(Cond, vecNmesh, "Nmesh", "z = "*string(z)*", Rmax = "*string(Rmax))
    (plt, Cond)
end

function conditionning(Rmax::Real, Nmesh::Int, vecz::AbstractArray, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64)
    Cond = zeros(T, length(vecz), 4)
    for (i, z) ∈ enumerate(vecz)
        Cond[i, :] = conditionning(Rmax, Nmesh, z, Basis, typemesh; opts_mesh, opts_basis, T)
    end
    plt = plot_conditionning(Cond, vecz, "z", "Nmesh = "*string(Nmesh)*", Rmax = "*string(Rmax))
    (plt, Cond)
end

function conditionning(Rmax::Real, Nmesh::Int, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64)
    # Default parameters
    Rmin = 0
    # Creation of the mesh
    m = typemesh(Rmin, Rmax, Nmesh; opts_mesh...)
    # Creation of the basis
    basis = Basis(m, T; opts_basis...)
    # Compute condition number
    conditionning(basis, z)
end

function conditionning(basis::Basis, z::Real)
    M₀  = mass_matrix(basis)
    A   = mass_matrix(deriv(basis))
    M₋₁ = weight_mass_matrix(basis, -1)
    H   = 1/2 * A - z .* M₋₁

    cond_M₀  = conditionning(M₀)
    cond_A   = conditionning(A)
    cond_M₋₁ = conditionning(M₋₁)
    cond_H   = conditionning(H)

    [cond_M₀, cond_A, cond_M₋₁, cond_H]
end

function conditionning(M::Matrix)
    pinvM = pinv(M)
    opnorm(M)*opnorm(pinvM)
end

function plot_conditionning(Cond::Matrix, X::AbstractVector, Xlabel::String, title::String)
    label = ["M₀", "A", "M₋₁", "H"]
    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :topright, 
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)     
    xlabel!(plt, Xlabel)
    ylabel!(plt, "Condition number") 
    title!(plt, title)
    for i ∈ 1:4
        plot!(plt, X, Cond[:,i], lw = 3, markershape = :x, markersize = 8, label = label[i], yaxis=:log)
    end
    plt
end