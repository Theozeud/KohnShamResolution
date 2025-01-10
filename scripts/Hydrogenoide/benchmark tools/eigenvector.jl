# This files provides function to compute the eigenvectors of the Hydrogenoid Schrodinger equations

function fundamental(z, x)
    z^(3/2)/sqrt(π)*exp(-z*abs(x))
end

function hydrogenoide_test_eigenvector(n::Int, vecRmax::AbstractVector, Nmesh::Int, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    eigs = []
    plts = Plots.Plot[]
    for Rmax ∈ vecRmax
        plt, eig = hydrogenoide_test_eigenvector(n, Rmax, Nmesh, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis, T = T, lₕ = lₕ)
        push!(eigs, eig)
        push!(plts, plt)
    end
    plot(plts..., layout = (length(plts)÷2, 2), size = (1200, 500 * length(plts)÷2))
end

function hydrogenoide_test_eigenvector(n::Int, Rmax::Real, vecNmesh::AbstractVector, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    eigs = []
    plts = Plots.Plot[]
    for Nmesh ∈ vecNmesh
        plt, eig = hydrogenoide_test_eigenvector(n, Rmax, Nmesh, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis, T = T, lₕ = lₕ)
        push!(eigs, eig)
        push!(plts, plt)
    end
    plot(plts..., layout = (length(plts)÷2, 2), size = (1200,500 * length(plts)÷2))
end

function hydrogenoide_test_eigenvector(n::Int, Rmax::Real, Nmesh::Int, vecz::AbstractVector, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    eigs = []
    plts = Plots.Plot[]
    for z ∈ vecz
        plt, eig = hydrogenoide_test_eigenvector(n, Rmax, Nmesh, z, Basis, typemesh; opts_mesh = opts_mesh, opts_basis = opts_basis, T = T, lₕ = lₕ)
        push!(eigs, eig)
        push!(plts, plt)
    end
    plot(plts..., layout = (length(plts)÷2, 2), size = (1200,500 * length(plts)÷2)), eigs
end

function hydrogenoide_test_eigenvector(n::Int, Rmax::Real, Nmesh::Int, z::Real, Basis, typemesh; opts_mesh = NamedTuple(), opts_basis = NamedTuple(), T = Float64, lₕ = 0)
    # Default parameters
    Rmin = 0
    # Creation of the mesh
    m = typemesh(Rmin, Rmax, Nmesh; opts_mesh...)
    # Creation of the basis
    basis = Basis(m, T; opts_basis...)
    # Creation of the Discretization
    D = KohnShamSphericalDiscretization(lₕ, basis, m)
    # Compute eigenvalues
    hydrogenoide_test_eigenvector(n, D, z, T)
end


function hydrogenoide_test_eigenvector(n::Int, D::KohnShamDiscretization, z::Real, T = Float64)

    method = ConstantODA(T(1))
    KM = KohnShamExtended(z = z, N = 1)

    sol = groundstate(KM, D, method; hartree = false, tol = 1e-16)

    fun = sol.eigvects[n]

    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :topright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
            
    xlabel!("r")
    ylabel!("Fundamental")
    title!("z = "*string(z)*", Nmesh = "*string(length(D.mesh))*", Rmax ="*string(last(D.mesh)))
    Rmin = 0
    X = range(Rmin, last(D.mesh)-0.0000001, 1000)
    
    if n == 1
        plot!(X, fundamental.(z,X),  label = "Théorique", color = :red, lw = 3)
    end
    
    plot!(plt, X, (sign(fun(X[1])) * fun).(X), label = "Numérique", color = :blue)

    (plt, sign(fun(X[1])) * fun)
end