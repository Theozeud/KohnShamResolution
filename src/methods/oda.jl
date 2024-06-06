struct ODA <: AbstractKohnShamResolutionMethod end

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

function init_cache(::ODA, model::AbstractDFTModel, discretization::KohnShamDiscretization)

    @unpack lₕ, Nₕ, basis = discretization

    # Init base matrices
    A   = zeros(Nₕ,Nₕ) # mass_matrix
    M₀  = zeros(Nₕ,Nₕ) # mass_matrix
    M₋₁ = zeros(Nₕ,Nₕ) # mass_matrix
    M₋₂ = zeros(Nₕ,Nₕ) # mass_matrix

    # Creation of the fix part of the hamiltonian   
    Kin =  zeros(lₕ+1, Nₕ, Nₕ)
    Coulomb =  zeros(lₕ+1, Nₕ, Nₕ)
    build_kinetic!(discretization, Kin, A, M₋₂)
    build_coulomb!(discretization, Coulomb, model, M₋₁)
    Hfix = Kin + Coulomb

    # Initialization of array for temporary stockage of computations
    temp_H       = zeros(Nₕ, Nₕ)

    temp_D       = zeros(lₕ+1, Nₕ, Nₕ)
    temp_Dstar   = zeros(lₕ+1, Nₕ, Nₕ)
    temp_U       = zeros(lₕ+1, Nₕ, Nₕ)
    temp_ϵ       = zeros(lₕ+1,Nₕ)
    temp_ϵ_sort  = zeros((lₕ+1)*Nₕ)
    temp_n       = zeros(lₕ+1, Nₕ)
    
    temp_tn      = 0.0     

    CacheODA(A, M₀, M₋₁, M₋₂, Hfix, temp_H, temp_Dstar, temp_D, temp_U, temp_ϵ, temp_ϵ_sort, temp_n, temp_tn)
end


function performstep!(::ODA, solver::KhonShamSolver)

    @unpack lₕ, Nₕ, basis = solver.discretization
    @unpack z, N, exc, potential = solver.model
    @unpack A, M₀, M₋₁, M₋₂, Hfix, temp_H, temp_Dstar, temp_D, temp_U, temp_ϵ, temp_ϵ_sort, temp_n, temp_tn = solver.cache

    # STEP 1 : find potential 
    #hartree = build_hartree!(discretization, matrix) 

    # STEP 2 : compute an approximation of the exchange correlation term
    #Exch = build_exchange_corr!(discretization, exchange_method, solver, matrix...)

    # STEP 3 : résolution du problème aux valeurs propres blocs par blocs
    for l ∈ 0:lₕ
        # Assembly the matrices
        temp_H .= Hfix[l+1] #+ Exch + Potential
        # Solve generalized eigenvalue problem on the section Hₗ
        temp_ϵ[l+1,:], temp_U[l+1,:,:] = solve_generalized_eigenvalue_problem(temp_H, M₀)
    end

    # STEP 4 : Build the n matrix using the Aufbau principle
    aufbau!(temp_n, temp_ϵ, N)

    # STEP 6 : Build the density related matrix
    build_density_star!(solver.discretization, temp_Dstar, temp_U, temp_n)

    # STEP 7 : update this matrix with a convex approach
    # some optimisation to find a good tₙ
    temp_tₙ = 0.5
    @. temp_D = temp_tₙ * temp_Dstar + (1 - temp_tₙ) * solver.Dprev

    # Registering into solver
    solver.U .= temp_U
    solver.n .= temp_n
    solver.ϵ .= temp_ϵ
    solver.D .= temp_D
end

stopping_criteria(m::ODA, solver::KhonShamSolver) = stopping_criteria(m, solver.D, solver.Dprev)
stopping_criteria(::ODA, D, Dprev) = norm(D .- Dprev)


function aufbau!(n::AbstractArray, ϵ::AbstractArray, N::Real)
    (d1,_) = size(ϵ)
    index_ϵ_sort = sortperm(ϵ, dims = 2)
    count = N
    i = firstindex(index_ϵ_sort)
    while count > 0
        l = div(index_ϵ_sort[i], d1)
        if count - (2*l + 1) ≥ 0
            for k in -l:l
                n[l+1,k+l+1] = 2
            end
            count = count - (2*l + 1)
        else
            for k in -l:l
                n[l+1,k+l+1] = 2/(2*l+1)
            end
            count = 0
            break
        end
        i+=1
    end
end


function update_density()


end

