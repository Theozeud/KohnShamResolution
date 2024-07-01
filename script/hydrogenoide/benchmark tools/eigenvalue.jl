m = logmesh(Rmin, Rmax, Nmesh; z = 0.5)
Rmax
basis = ShortP1Basis(m, T;  normalize = true, left = false, right = false)

lₕ = 0
z




function hydrogenoide_test_eigenvalue(basis::Basis, mesh::OneDMesh, z::Real, lₕ::Int, Rmax::Real)

    T = bottom_type(basis)

    function eigval_theo(n, z)
        -z^2/(2*n^2)
    end
    
    D = KohnShamSphericalDiscretization(lₕ, basis, mesh)
    KM = KohnShamExtended(z = z, N = 1)

    method = ConstantODA(T(1))
    sol = groundstate(KM, D, method; hartree = false)

    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :bottomright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)

    xlabel!("n")
    ylabel!("Energie")
    title!("z = "*string(z)*", Nmesh = "*length(mesh)*", Rmax ="*string(Rmax))

    index_ϵ = findall(x->x < 0, sol.ϵ)

    scatter!(index_ϵ, eigval_theo.(index_ϵ,z),  label = "Théorique",
                markershape = :circ, 
                markersize = 8,
                markeralpha = nothing,
                markercolor = :red,
                markerstrokewidth = 1,
                markerstrokecolor = :red
                )

    scatter!(index_ϵ, sol.ϵ[index_ϵ],  label = "Numérique",
                markershape = :x, 
                markersize = 10,
                markeralpha = nothing,
                markercolor = :black,
                markerstrokewidth = 2)

    plt
end