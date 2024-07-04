# This files provides function to compute the order of convergence of methods applied to the Hydrogenoid Schrodinger equations

function eigval_theo(n, z)
    -z^2/(2*n^2)
end

function fundamental(z, x)
    z^(3/2)/sqrt(π)*exp(-z*abs(x))
end

#Rmax = (1.5 * log(z) + cutting_pre*log(10))/z

function test_convergence_withNmesh(vecNmesh::AbstractVector, Rmax::Real, z::Real, vecBasis::NamedTuple, typemesh; opts_mesh = NamedTuple(), opts_basis = [NamedTuple() for i ∈ eachindex(vecBasis)], T = Float64, lₕ = 0, nb_eigval = 1)
    plt_ϵ_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the first "*string(nb_eigval)*"-th eigenvalues")
    title!(plt_ϵ_error, "Rmax = "*string(Rmax)*" and z = "*string(z))

    plt_u_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_u_error, "Nmesh")
    ylabel!(plt_u_error, "L2 error on the fundamental")
    title!(plt_u_error, "Rmax = "*string(Rmax)*" and z = "*string(z))

    Uerror = zeros(length(vecBasis), length(vecNmesh))
    ϵerror = zeros(length(vecBasis), length(vecNmesh), nb_eigval)
    for (i,Basis) ∈ enumerate(vecBasis)
        println("Basis "*string(keys(vecBasis)[i]))
        (_, _, ϵ_error, u_error) = test_convergence_withNmesh(vecNmesh, Rmax, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis[i], T = T, lₕ = lₕ, nb_eigval = nb_eigval)
        for j ∈ 1:nb_eigval
            plot!(plt_ϵ_error, vecNmesh, ϵ_error[:, j], lw = 3, label = "Basis "*string(keys(vecBasis)[i])*", ϵ"*string(j), markershape = :x, markersize = 10)
        end
        plot!(plt_u_error, vecNmesh, u_error, lw = 3, label = "Basis "*string(keys(vecBasis)[i]), markershape = :x, markersize = 10)
        Uerror[i, :] = u_error
        ϵerror[i, : , :] = ϵ_error
    end
    (plt_ϵ_error, plt_u_error, ϵerror, Uerror)
end


function test_convergence_withNmesh(vecNmesh::AbstractVector, Rmax::Real, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0, nb_eigval = 1)
    # Default parameters
    Rmin = zero(typeof(Rmax))
    # Setup the model
    KM = KohnShamExtended(z = z, N = 1)
    # Setup the method
    method = ConstantODA(one(T))
    # Computation of the errors
    ϵ_error = zeros(T, length(vecNmesh), nb_eigval)
    u_error = zeros(T, length(vecNmesh))
    for (i,Nmesh) ∈ enumerate(vecNmesh)
        # Creation of the discretization
        m = typemesh(Rmin, Rmax, Nmesh; opts_mesh...)
        #@show opts_basis
        basis = Basis(m, T; opts_basis...)
        D = KohnShamSphericalDiscretization(lₕ, basis, m)
        # Solving the problem
        @time "Nmesh = "*string(Nmesh) sol = groundstate(KM, D, method; tol = 1e-20, hartree = false)
        funda = sol.eigvects[1]
        # Compute the error for eigenvalues and the fundamental
        ϵ_error[i,:] = abs.(sol.ϵ[1:nb_eigval] .- eigval_theo.(1:nb_eigval, z))
        u_error[i]   = sqrt(1/(Nmesh - 1) * sum((funda.(m.points[1:end-1]) .- fundamental.(z, m.points[1:end-1])).^2))
    end
    # Creation of the plot for eigenvalue
    plt_ϵ_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
                                 legendfontsize  = 12,  
                                 titlefontsize   = 12,
                                 guidefontsize   = 12,
                                 tickfontsize    = 12)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the first "*string(nb_eigval)*"-th eigenvalues")
    title!(plt_ϵ_error, "Rmax = "*string(Rmax)*" and z = "*string(z))
    for i ∈ 1:nb_eigval
        plot!(plt_ϵ_error, vecNmesh, ϵ_error[:, i], lw = 3, label = "ϵ"*string(i), markershape = :x, markersize = 10)
    end
    # Creation of the plot for the fundamental
    plt_u_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
                                 legendfontsize  = 12,  
                                 titlefontsize   = 12,
                                 guidefontsize   = 12,
                                 tickfontsize    = 12)
    xlabel!(plt_u_error, "Nmesh")
    ylabel!(plt_u_error, "L2 error on the fundamental")
    title!(plt_u_error, "Rmax = "*string(Rmax)*" and z = "*string(z))
    plot!(plt_u_error, vecNmesh, u_error, lw = 3, label = "Fundamental ", markershape = :x, markersize = 10)
    return (plt_ϵ_error, plt_u_error, ϵ_error, u_error)
end
