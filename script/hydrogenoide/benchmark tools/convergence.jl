# This files provides function to compute the order of convergence of methods applied to the Hydrogenoid Schrodinger equations

function eigval_theo(n, z)
    -z^2/(2*n^2)
end

function fundamental(z, x)
    z^(3/2)/sqrt(π)*exp(-z*abs(x))
end

function test_convergence_withNmesh(vecNmesh::AbstractVector, Rmax::Real, z::Real, vecBasis::NamedTuple, typemesh; opts_mesh = NamedTuple(), opts_basis = [NamedTuple() for i ∈ eachindex(vecBasis)], T = Float64, l = 0, nb_eigval = 1)
    plt_ϵ_error = plot( size = (750,550), margin = 0.5Plots.cm, legend = :outertopright, xaxis=:log, yaxis=:log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the $nb_eigval-th eigenvalue")
    title!(plt_ϵ_error, "Rmax = $Rmax, z = $z and l = $l")

    ϵerror = zeros(length(vecBasis), length(vecNmesh))
    for (i,Basis) ∈ enumerate(vecBasis)
        println("Basis "*string(keys(vecBasis)[i]))
        (_, ϵ_error) = test_convergence_withNmesh(vecNmesh, Rmax, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis[i], T = T, l = l, nb_eigval = nb_eigval)
        plot!(plt_ϵ_error, vecNmesh, ϵ_error, lw = 4, label = string(keys(vecBasis)[i]), markershape = :x, markersize = 10)
        ϵerror[i, : ] = ϵ_error
    end
    (plt_ϵ_error, ϵerror)
end


function test_convergence_withNmesh(vecNmesh::AbstractVector, Rmax::Real, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, l = 0, nb_eigval = 1)
    # Default parameters
    Rmin = zero(typeof(Rmax))
    # Setup the model
    KM = KohnShamExtended(z = z, N = 1)
    # Setup the method
    method = ConstantODA(one(T))
    # Computation of the errors
    ϵ_error = zeros(T, length(vecNmesh))
    for (i,Nmesh) ∈ enumerate(vecNmesh)
        # Creation of the discretization
        m = typemesh(Rmin, Rmax, Nmesh; opts_mesh...)
        basis = Basis(m, T; opts_basis...)
        D = KohnShamRadialDiscretization(l, basis, m)
        # Solving the problem
        @time "Nmesh = $Nmesh" sol = groundstate(KM, D, method; tol = 1e-20, hartree = false)
        # Compute the error for eigenvalues
        ϵ_error[i] = abs(sol.ϵ[nb_eigval] - eigval_theo(nb_eigval+l, z))
    end
    # Creation of the plot for eigenvalue
    plt_ϵ_error = plot( size = (1300,1000), margin = 0.5Plots.cm, legend = :outertopright, xaxis=:log, yaxis=:log,
                                 legendfontsize  = 12,  
                                 titlefontsize   = 12,
                                 guidefontsize   = 12,
                                 tickfontsize    = 12)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the $nb_eigval-th eigenvalues")
    title!(plt_ϵ_error, "Rmax = $Rmax, z = $z and l = $l")
    plot!(plt_ϵ_error, vecNmesh, ϵ_error, lw = 4, label = "ϵ$nb_eigval", markershape = :x, markersize = 10)
    return (plt_ϵ_error, ϵ_error)
end
