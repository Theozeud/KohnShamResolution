using KohnShamResolution
using LinearAlgebra
using Plots

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0
Rmin = 0

# One electron model
z = 1
N = 1
KM = KohnShamExtended(z = z, N = N) 
eigvalue_theo(n,z)    = -z^2/(2*n^2)
fundamental(z, x)     = z^(3/2)/sqrt(π)*exp(-z*abs(x))

# Vector of Error
err_eigenvector = []
err_eigenvalue = []
condition_number = []
RmaxArray = 10:2:150
#cutting_preArray = [5,10,15]

for Rmax in RmaxArray

    #Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
    #push!(RmaxArray, Rmax)
    m = logmesh(0,Rmax,100)
    basis = HatBasis(m; left = false, right = false)
    D = KohnShamSphericalDiscretization(lₕ, basis, m)

    # Computation by hand
    deriv_basis = deriv(basis)

    A   = mass_matrix(deriv_basis, Rmin, Rmax)
    M₀  = mass_matrix(basis, Rmin, Rmax)
    M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
    M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

    H = 1/2 * (A + lₕ*(lₕ+1)*M₋₂) - z .* M₋₁

    ϵ, U = eigen(H,M₀)

    eig1 = KohnShamResolution.build_on_basis(basis, U[1,:]) * Monomial(-1)
    fun =  eig1 / sqrt(scalar_product(eig1,eig1, m))  

    push!(err_eigenvector, norm() / 98)
    push!(err_eigenvalue, 100 * abs(ϵ[begin] - eigvalue_theo(1,z))/abs(eigvalue_theo(1,z)) )
    push!(condition_number, abs(ϵ[end]/ϵ[begin]))

end


# Plots

plt_eigenvector = plot( RmaxArray, err_eigenvector,  size = (900,600), margin = 0.5Plots.cm, legend = :bottomright, 
                        lw = 3, color = :black, markershape = :x, markersize = 8,
                        legendfontsize  = 14,  
                        titlefontsize   = 14,
                        guidefontsize   = 14,
                        tickfontsize    = 14)

xlabel!(plt_eigenvector, "Rmax")
ylabel!(plt_eigenvector, "Error on the fundamental")
title!(plt_eigenvector, "z = "*string(z))

plt_eigenvalue = plot(  RmaxArray, err_eigenvalue, size = (900,600), margin = 0.5Plots.cm, legend = :bottomright, 
                        lw = 3, color = :black, markershape = :x, markersize = 8,
                        legendfontsize  = 14,  
                        titlefontsize   = 14,
                        guidefontsize   = 14,
                        tickfontsize    = 14)

xlabel!(plt_eigenvalue, "Rmax")
ylabel!(plt_eigenvalue, "Error on the minimal eigenvalue (%)")
title!(plt_eigenvalue, "z = "*string(z))

plt_cn = plot(  RmaxArray, condition_number, size = (900,600), margin = 0.5Plots.cm, legend = :bottomright, 
                lw = 3, color = :black, markershape = :x, markersize = 8,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
                
xlabel!(plt_cn, "Rmax")
ylabel!(plt_cn, "Conditionning number")
title!(plt_cn, "z = "*string(z))


pltfin = plot(plt_eigenvector, plt_eigenvalue, plt_cn, layout = (1,3), size = (2400,1000), margin = 1.2Plots.cm)
savefig(pltfin, "image/hydrogenoide/with P1 elements/Study of the impact of Rmax")