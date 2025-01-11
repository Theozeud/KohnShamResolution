
#############################################################################################
#                                     PLOT EIGENVECTOR
#############################################################################################


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