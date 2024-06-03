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
    temp_n
    temp_ϵ
    temp_U
    temp_Rstar
    temp_R
    temp_tn
end

function init_cache(model::AbstractDFTModel,::ODA; lₕ, Nₕ)

    # Init base matrices
    A   = zeros(Nₕ,Nₕ)
    M₀  = zeros(Nₕ,Nₕ)
    M₋₁ = zeros(Nₕ,Nₕ)
    M₋₂ = zeros(Nₕ,Nₕ)

    # Creation of the fix part of the hamiltonian for each section Hₗ    
    Hfix = zeros(lₕ+1, Nₕ, Nₕ) 
    for l ∈ 0:lₕ
        Hfix[l+1,:,:] .= - A .- l*(l+1)*M₋₂ .- 2*charge(model)*(2*l+1) .* M₋₁
    end 

    # Initialization of array for temporary stockage of computations
    temp_H       = zeros(Nₕ, Nₕ)
    temp_U       = zeros(lₕ+1, (2lₕ+1)Nₕ, Nₕ)
    temp_n       = zeros(lₕ+1, (2lₕ+1)Nₕ)
    temp_ϵ       = zeros(lₕ+1,(2lₕ+1)Nₕ)
    temp_R       = zeros(lₕ+1, Nₕ, Nₕ)
    temp_Rstar   = zeros(lₕ+1, Nₕ, Nₕ)
    temp_tn      = 0.0     

    CacheODA(A, M₀, M₋₁, M₋₂, Hfix, temp_H, temp_n, temp_ϵ, temp_U, temp_Rstar, temp_R, temp_tn)
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
    for l ∈ 0:lₕ   
        for k ∈ 1:(2*l+1)*Nₕ
            if cache.temp_ϵ[l,k] != ϵf
                temp_n[l,k] = cache.temp_ϵ[l,k] < ϵf ? 2 : 0
            else
                
            end
        end
    end

    # STEP 6 : Build the density related matrix 
    cache.temp_Rstar = zeros(lₕ+1,Nₕ,Nₕ)
    for l ∈ 0:lₕ
        for k ∈ 1:(2*l+1)*Nₕ
            cache.temp_Rstar[l] .+= cache.temp_n[l,k]*tensorproduct(cache.temp_U[l,k], cache.temp_U[l,k])
        end
    end

    # STEP 7 : update this matrix with a convex approach
    # some optimisation to find a good tₙ
    cache.temp_R .= cache.temp_tₙ .* cache.temp_Rstar + (1 - cache.temp_tₙ) .* solver.Rprev

    # Registering into solver
    solver.U .= cache.temp_U
    solver.n .= cache.temp_n
    solver.ϵ .= cache.temp_ϵ
    solver.R .= cache.temp_R
end

function stopping_criteria(::ODA, R, Rprev)
    val = abs(R .- Rprev)
    solver.current_crit = val
    val
end

stopping_criteria(m::ODA, solver::KhonShamSolver) = stopping_criteria(m, solver.R, solver.Rprev)