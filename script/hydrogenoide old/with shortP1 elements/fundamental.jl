using KohnShamResolution
using Plots

# Choice of the method
method = ConstantODA(1.0)

# Discretization 
lₕ = 0
Nmesh = 100
Rmin = 0.00001

# One electron model
zA = [1, 2, 3, 4]
N = 1
eigvalue_theo(n,z)    = -z^2/(2*n^2)
fundamental(z, x)     = z^(3/2)/sqrt(π)*exp(-z*abs(x))

# Plot
pltA = []

cutting_pre = 5

for z in zA

    Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
    m = logmesh(Rmin, Rmax, Nmesh, z = 0.5)
    basis = ShortP1Basis(m; left = false, right = false, normalize = true)
    D = KohnShamSphericalDiscretization(lₕ, basis, m)

    KM = KohnShamExtended(z = z, N = N)

    sol = groundstate(KM, D, method; tol = 1e-20, hartree = false, maxiter = 1)

    fun = sol.eigvects[1]

    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :bottomright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)

    xlabel!("r")
    ylabel!("Fundamental")
    title!("z = "*string(z))

    X = logmesh(Rmin, Rmax,1000).points

    plot!(X, fundamental.(z,X),  label = "Théorique", color = :red, lw = 3)

    plot!(plt, X, -fun.(X), label = "Numérique", color = :blue)

    push!(pltA, plt)
end

pltfin = plot(pltA..., layout = (2,2), size = (1200,1000))
savefig(pltfin, "image/hydrogenoide/with shortP1 elements/Vecteurs propres avec Nmesh ="*string(Nmesh))