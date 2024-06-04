struct ODA <: AbstractKohnShamResolutionMethod end

for_sphericalsymmetry(::ODA) = true
for_cylindricalsymmetry(::ODA) = false
ismethod_for_model(::ODA,::KohnShamExtended) = true

struct CacheODA <: AbstractKohnShamCache
    A
    M₀
    M₋₁
    M₋₂
    Hfix
    temp_H
    temp_Dstar
    temp_D
    temp_U
    temp_ϵ
    temp_ϵ_sort
    temp_n    
    temp_tn
end

function init_cache(::ODA, model::AbstractDFTModel, discretization::KohnShamDiscretization, ; lₕ, Nₕ)

    @unpack lₕ, Nₕ, basis = discretization

    # Init base matrices
    A   = zeros(Nₕ,Nₕ) # mass_matrix
    M₀  = zeros(Nₕ,Nₕ) # mass_matrix
    M₋₁ = zeros(Nₕ,Nₕ) # mass_matrix
    M₋₂ = zeros(Nₕ,Nₕ) # mass_matrix

    # Creation of the fix part of the hamiltonian for each section Hₗ    
    Hfix = zeros(lₕ+1, Nₕ, Nₕ) 
    for l ∈ 0:lₕ
        Hfix[l+1,:,:] .= - A .- l*(l+1)*M₋₂ .- 2*charge(model)*(2*l+1) .* M₋₁
    end 

    # Initialization of array for temporary stockage of computations
    temp_H       = zeros(Nₕ, Nₕ)

    temp_D       = zeros(lₕ+1, Nₕ, Nₕ)
    temp_Dstar   = zeros(lₕ+1, Nₕ, Nₕ)
    temp_U       = zeros(lₕ+1, (2lₕ+1)*Nₕ, Nₕ)
    temp_ϵ       = zeros(lₕ+1,(2lₕ+1)*Nₕ)
    temp_ϵ_sort  = zeros((lₕ+1)*(2lₕ+1)*Nₕ)
    temp_n       = zeros(lₕ+1, (2lₕ+1)*Nₕ)
    
    temp_tn      = 0.0     

    CacheODA(A, M₀, M₋₁, M₋₂, Hfix, temp_H, temp_Dstar, temp_D, temp_U, temp_ϵ, temp_ϵ_sort, temp_n, temp_tn)
end


function performstep!(::ODA, cache::CacheODA , solver::KhonShamSolver)

    # STEP 1 : find potential 
    #Potential = 

    # STEP 2 : compute an approximation of the exchange correlation term
    #Exch = 

    # STEP 3 : résolution du problème aux valeurs propres blocs par blocs
    for l ∈ 0:lₕ
        # Assembly the matrices
        cache.temp_H .= cache.Hix[l+1] #+ Exch + Potential
        # Solve generalized eigenvalue problem on the section Hₗ
        cache.temp_ϵ[l+1], cache.temp_U[l] = solve_generalized_eigenvalue_problem(cache.temp_H, cache.M₀)
    end

    # STEP 4 : Build the n matrix using the Aufbau principle
    aufbau!(cache.temp_n, cache.temp_ϵ, z)

    # STEP 6 : Build the density related matrix 
    cache.temp_Rstar = zeros(lₕ+1,Nₕ,Nₕ)
    for l ∈ 0:lₕ
        for k ∈ 1:(2*l+1)*Nₕ
            cache.temp_Rstar[l] .+= cache.temp_n[l,k]*tensorproduct(cache.temp_U[l,k], cache.temp_U[l,k])
        end
    end

    # STEP 7 : update this matrix with a convex approach
    # some optimisation to find a good tₙ
    cache.temp_tₙ = 0.5
    cache.temp_R .= cache.temp_tₙ .* cache.temp_Rstar + (1 - cache.temp_tₙ) .* solver.Rprev

    # Registering into solver
    solver.U .= cache.temp_U
    solver.n .= cache.temp_n
    solver.ϵ .= cache.temp_ϵ
    solver.R .= cache.temp_R
end

function stopping_criteria(::ODA, R, Rprev)
    norm(R .- Rprev)
end

stopping_criteria(m::ODA, solver::KhonShamSolver) = stopping_criteria(m, solver.R, solver.Rprev)


function aufbau!(n::AbstractArray, ϵ::AbstractArray, N::Real)
    index_ϵ_sort = sortperm(ϵ)
    count = N
    i = firstindex(index_ϵ_sort)
    while count > 0
        l = div(index_ϵ_sort[i], Nₕ)
        if count - (2*l + 1) ≥ 0
            for k in -l:l
                n[l,k] = 2
            end
            count = count - (2*l + 1)
        else
            for k in -l:l
                n[l,k] = 2/(2*l+1)
            end
            count = 0
            break
        end
        i+=1
    end
end


