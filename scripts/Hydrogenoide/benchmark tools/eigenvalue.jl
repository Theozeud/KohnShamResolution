# This files provides function to compute the eigenvalues of the Hydrogenoid Schrodinger equations

function eigval_theo(n, z)
    -z^2/(2*n^2)
end

function hydrogenoide_test_eigenvalue(vecRmax::AbstractVector, Nmesh::Int, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    eigs = Vector{T}[]
    plts = Plots.Plot[]
    for Rmax ∈ vecRmax
        plt, eig = hydrogenoide_test_eigenvalue(Rmax, Nmesh, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis, T = T, lₕ = lₕ)
        push!(eigs, eig)
        push!(plts, plt)
    end
    plot(plts, layout = (length(plts)÷2+1, 2), size = (1200,1000))
end

function hydrogenoide_test_eigenvalue(Rmax::Real, vecNmesh::AbstractVector, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    eigs = Vector{T}[]
    plts = Plots.Plot[]
    for Nmesh ∈ vecNmesh
        plt, eig = hydrogenoide_test_eigenvalue(Rmax, Nmesh, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis, T = T, lₕ = lₕ)
        push!(eigs, eig)
        push!(plts, plt)
    end
    plot(plts, layout = (length(plts)÷2+1, 2), size = (1200,1000))
end

function hydrogenoide_test_eigenvalue(Rmax::Real, Nmesh::Int, vecz::AbstractVector, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    eigs = Vector{T}[]
    plts = Plots.Plot[]
    for z ∈ vecz
        plt, eig = hydrogenoide_test_eigenvalue(Rmax, Nmesh, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis, T = T, lₕ = lₕ)
        push!(eigs, eig)
        push!(plts, plt)
    end
    plot(plts..., layout = (length(plts)÷2, 2), size = (1200,1000)), eigs
end

function hydrogenoide_test_eigenvalue(Rmax::Real, Nmesh::Int, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    # Default parameters
    Rmin = 0
    # Creation of the mesh
    m = typemesh(Rmin, Rmax, Nmesh; opts_mesh...)
    # Creation of the basis
    basis = Basis(m, T; opts_basis...)
    # Creation of the Discretization
    D = KohnShamSphericalDiscretization(lₕ, basis, m)
    # Compute eigenvalues
    hydrogenoide_test_eigenvalue(D, z, T)
end


function hydrogenoide_test_eigenvalue(D::KohnShamDiscretization, z::Real, T = Float64)

    method = ConstantODA(T(1))
    KM = KohnShamExtended(z = z, N = 1)

    sol = groundstate(KM, D, method; hartree = false, tol = 1e-16)

    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :bottomright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)

    xlabel!("n")
    ylabel!("Energie")
    title!("z = "*string(z)*", Nmesh = "*string(length(D.mesh))*", Rmax ="*string(Rmax))

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

    (plt, sol.ϵ[index_ϵ])
end