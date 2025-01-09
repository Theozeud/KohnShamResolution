# Plot convergence avec Nmesh et Rmax !
    # en terme de valeurs propres
    # en terme de vecteurs propres
        # -> il faut une norme


# Courbe convergence avec Nmesh et Rmax
# On veut pouvoir comparer pour différentes bases et mesh (typebase, typemesh)

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



# z,L

# basis <- mesh <- Rmax, Nmesh, typemesh, opts_mesh
#       <- T
#       <- opts_basis