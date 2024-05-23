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
    temp_U
    temp_Rstar
end

function init_cache(::ODA)

    # Init base matrices
    #A
    #M₀
    #M₋₁
    #M₋₂

    # Creation of the Hamiltonian for each section Hₗ    
    H = zeros(lₕ, Nₕ, Nₕ) 
    for l ∈ 0:lₕ
        H[l+1,:,:] .= - A .- l*(l+1)*M₋₂ .- 2*z*(2*l+1).*M₋₁
    end




end


function performstep!(method::ODA, cache::CacheODA , solver::KhonShamSolver)

    # STEP 1 : find potential 
    #Potential = 

    # STEP 2 : compute an approximation of the exchange correlation term
    #Exch = 

    # STEP 3 : résolution du problème aux valeurs propres blocs par blocs
    for l ∈ 0:lₕ
        # Assembly the matrices
        cache.temp_H .= cache.H[l+1] #+ Exch + Potential
        # Solve generalized eigenvalue problem on the section Hₗ
        cache.temp_U, cache.temp_ϵ .= solve_generalized_eigenvalue_problem(temp_H, cache.M₀)
    end

    # STEP 4 : Building the n matrix using the Aufbau principle
    for l ∈ 0:lₕ   
        for k ∈ 1:(2*l+1)*Nₕ
            if cache.temp_ϵ[l,k] != ϵf
                temp_n[l,k] = cache.temp_ϵ[l,k] < ϵf ? 2 : 0
            else

            end
        end
    end

    # STEP 5 : Build the density related matrix 
    cache.temp_Rstar = zeros(lₕ+1,Nₕ,Nₕ)
    for l ∈ 0:lₕ
        for k ∈ 1:(2*l+1)*Nₕ
            cache.temp_Rstar[l] .+= cache.temp_n[l,k]*tensorproduct(cache.temp_U[l,k], cache.temp_U[l,k])
        end
    end

    # STEP 6 : update this matrix with a convex approach



    # Registering into solver
    solver.U .= U
    solver.n .= n
    solver.ϵ .= ϵ
    solver.R .= R
end

function stopping_criteria(method, solver)


end