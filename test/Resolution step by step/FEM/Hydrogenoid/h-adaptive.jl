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
    Nm = length(m)
    error = zeros(T, Nm - 1)
    eval = eval_sol2(U[:,1], basis)
    for i ∈ 1:Nm-1
        error[i] = abs(eval(m[i]) - eigvect_theo(T, 1, z)(m[i]))/sqrt((m[i+1] - m[i]))
    end
    error
end

# h-strategy
function h_strategy(m, localerr, tol, Nmax)
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


# full procedure
function hadapt_hydro(mesh, ordermax, tol, Nmax; l = 0, T = Float64, z = 1)
    Mesh = []
    refine = true
    iter = 1
    Λ, U = zero(T), zeros(1,1)
    basis = nothing
    while refine
        println("iter = $iter")
        @show mesh
        push!(Mesh, mesh)
        nU, Λ, basis =  eigen_hydro(mesh, ordermax, l, T)
        localerr = local_errors(Λ, nU, basis, mesh, z, T)
        mesh, refine = h_strategy(mesh, localerr, tol, Nmax)
        iter += 1
        refine = false
        U = nU
    end
    return Λ, U, Mesh, basis
end

l = 0
z = 1
Rmax = 10
Nmbegin = 10
tol = 1e-8
Nmax = 80
order = 3


@show m = mesh(Float64.([0,0.5,1,2,3,4,5,6,7,8,9,10]))
#linmesh(0,Rmax,Nmbegin)

sol = hadapt_hydro(m, order, tol, Nmax, l = l, z = z)
Λ = sol[1]
U = sol[2]
Mesh = sol[3]
basis = sol[4]

@show abs(-0.5 -Λ[1])

using Plots
unum2 = eval_sol2(U[:,1], basis)
X = LinRange(0,Rmax, 100)
plot(X, eigvect_theo(Float64, 1, z).(X), lw = 3, label = "Théorique")
plot!(X, abs.(unum2.(X)), label = "h-adaptative", lw = 3)
plot!(X, abs.([eval_sol(U, basis, x) for x ∈ X]), lw = 3, markershape = :x)


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


