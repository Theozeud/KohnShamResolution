using KohnShamResolution
using LinearAlgebra
using Plots


# theoretical Eigenvalue and Eigenvectors
function eigval_theo(T, n, z)
    -T(z)^2/(T(2)*T(n)^2)
end

function eigvect_theo(T, n, z)
    if n == 1
        return x->T(z)^(3/2)/sqrt(π)*exp(-z*abs(x))*x
    end
end

# eigenvalues and eigenvectors with IntLeg
function L2normalization!(U, M₀)
    ax = axes(U, 1)
    nU = zero(U)
    for k ∈ axes(U,2)
        normalization = sqrt(sum([U[i,k] * U[j,k] * M₀[i,j] for i∈ax for j∈ax]) * 4π)
        nU[:,k] .= U[:,k] / normalization
    end
    nU
end


function eigen_hydro(mesh, ordermax, l = 0, T = Float64)
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = 2, ordermax = ordermax,  normalize = false)
    deriv_basis = deriv(basis)
    A   = Symmetric(mass_matrix(deriv_basis))
    M₀  = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    if l == 0
        H = T(0.5) * A  -  M₋₁
        Λ, U = eigen(H, M₀)
        nU = L2normalization!(U, M₀)
        return  (nU, Λ, basis)
    else
        M₋₂ = Symmetric(weight_mass_matrix(basis, -2))
        H = T(0.5) * A  -  M₋₁ + T(0.5) * l*(l+1) * M₋₂
        Λ, U = eigen(H, M₀)
        nU = L2normalization!(U, M₀)
        return  (nU, Λ, basis)
    end
end

# computation of local errors
function eval_sol(U, basis, x)
    ux = 0
    for i ∈ axes(U,2)
        ux += U[i] * KohnShamResolution.eval_basis(basis, i, x)
    end
    ux
end

function eval_sol2(U, basis)
    build_on_basis(basis, U)
end

function local_errors(Λ, U, basis, m, z = 1, T = Float64)
    Nm = length(m.points)
    error = zeros(T, Nm - 1)
    eval = eval_sol2(U[:,1], basis)
    for i ∈ 1:Nm-1
        error[i] = abs(eval(m[i+1]) - eigvect_theo(T, 1, z)(m[i+1]))/sqrt((m[i+1] - m[i])) + abs(Λ[1] - eigval_theo(T, 1, z))
    end
    error
end

# h-strategy
function h_strategy1(m, localerr, tol, Nmax)
    points = copy(m.points)
    Nm = length(m)
    index = sortperm(localerr; rev = true)
    cumul2err = zeros(eltype(localerr), Nm-2)
    for i ∈ 2:Nm-1
        cumul2err[i-1] = sqrt(localerr[i-1]^2+localerr[i]^2)
    end
    index2 = sortperm(cumul2err; rev = true)
    newpoints = [first(points), last(points)]
    Nb = 2
    while Nb < Nmax
        if !isempty(index) && !isempty(index2)
            fi = first(index)
            fi2 = first(index2)
            if localerr[fi] > cumul2err[fi2] && localerr[fi] > tol/Nm
                push!(newpoints, (points[fi] + points[fi+1])/2)
                popfirst!(index)
            elseif localerr[fi] ≤ cumul2err[fi2] && cumul2err[fi2] > tol/Nm
                push!(newpoints, points[fi2+1])
                popfirst!(index2)
            else
                break
            end
        elseif !isempty(index)
            fi = first(index)
            if  localerr[fi] > tol/Nmax
                push!(newpoints, (points[fi] + points[fi+1])/2)
                popfirst!(index)
            else
                break
            end
        elseif !isempty(index2)
            fi2 = first(index2)
            if  cumul2err[fi2] > tol/Nmax
                push!(newpoints, points[fi2+1])
                popfirst!(index2)
            else
                break
            end
        else
            break
        end
        Nb += 1
    end
    mesh(newpoints), (Nb < Nmax) 
end


function h_strategy2(m, localerr, tol, Nmax)
    newpoints = copy(m.points)
    Nm = length(m)
    index = sortperm(localerr; rev = true)
    Nb = Nm
    while Nb < Nmax && !isempty(index)
        fi = first(index)
        if localerr[fi] > tol/Nm
            push!(newpoints, (m.points[fi] + m.points[fi+1])/2)
            popfirst!(index)
            Nb += 1
        else
            break
        end
    end
    mesh(newpoints), (Nb < Nmax) 
end


# full procedure
function hadapt_hydro(m, ordermax, tol, Nmax; l = 0, T = Float64, z = 1)
    APlot = [plot_sol_inter(0, m, m)]
    Mesh = [m]
    refine = true
    iter = 1
    Λ, U = zero(T), zeros(1,1)
    basis = nothing
    while refine
        println("iter = $iter")
        nU, Λ, basis =  eigen_hydro(m, ordermax, l, T)
        localerr = local_errors(Λ, nU, basis, m, z, T)
        m, refine = h_strategy2(m, localerr, tol, Nmax)
        plt = plot_sol_inter(iter, last(Mesh), m)
        push!(APlot, plt)    
        push!(Mesh, m)
        iter += 1
        U = nU
        refine = false
    end
    return Λ, U, Mesh, basis, APlot
end


# Plot intermédiare
function plot_sol_inter(step::Int, last_mesh, current_mesh)
    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :outertopright,
                legendfontsize  = 18,  
                titlefontsize   = 18,
                guidefontsize   = 18,
                tickfontsize    = 18)
    xlabel!(plt, "t")
    title!(plt, "Step $step")
    scatter!(plt, current_mesh.points, zero(current_mesh.points); markershape = :o, markercolor = :red, markersize = 12)
    scatter!(plt, last_mesh.points, zero(last_mesh.points); markershape = :o, markercolor = :black, markersize = 12)
    plt
end



## Script

l = 1
z = 1
Rmax = 80
Nmbegin = 90
tol = 1e-4
Nmax = 20
order = 4


@show m = geometricmesh(0, Rmax, Nmbegin;s = 0.9) #mesh(Float64.([0,1,2,3,4,5,6,7,8,9,10]))
#linmesh(0, Rmax, Nmbegin) #

sol = hadapt_hydro(m, order, tol, Nmax, l = l, z = z)
Λ = sol[1]
U = sol[2]
Mesh = sol[3]
basis = sol[4]
APlot = sol[5]

@show abs(-0.125 -Λ[1])

using Plots
unum2 = eval_sol2(U[:,1], basis)
X = LinRange(0,Rmax, 1000)
plot(size = (1000,700))
plot!(X, eigvect_theo(Float64, 1, z).(X), lw = 3, label = "Théorique")
plot!(X, abs.(unum2.(X)), label = "h-adaptative", lw = 3)
#plot!(X, abs.([eval_sol(U, basis, x) for x ∈ X]), lw = 3, markershape = :x)


#=
@show m = mesh(Float64.([0,1,2,4,6,8,10]))
basis = ShortP1IntLegendreBasis(m, Float64; ordermin = 2, ordermax = order,  normalize = false)

################
# Plots of basis
X = LinRange(0, Rmax, 10000)
plt_basis = plot(legend = false)
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(plt_basis, X, p.(X), lw = 3)
end
xlabel!("Maillage")
plt_basis

#######################################
# Plots of the derivatives of the basis
X = LinRange(0, Rmax, 10000)
deriv_basis = deriv(basis)
plt_derivbasis = plot(legend = false)
for i ∈ 1:length(deriv_basis)
    p = build_basis(deriv_basis, i)
    plot!(plt_derivbasis, X, p.(X), lw = 3)
end
plt_derivbasis


A   = Symmetric(mass_matrix(deriv_basis))
M₀  = Symmetric(mass_matrix(basis))
=#


