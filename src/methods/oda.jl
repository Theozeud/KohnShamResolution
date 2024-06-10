struct ODA <: AbstractKohnShamResolutionMethod end

struct CacheODA <: AbstractKohnShamCache
    A
    M₀
    M₋₁
    M₋₂
    Hfix
    tmp_H
    tmp_Dstar
    tmp_D
    tmp_U
    tmp_ϵ
    tmp_ϵ_sort
    tmp_n    
    tmp_tn
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
    tmp_H       = zeros(Nₕ, Nₕ)

    tmp_D       = zeros(lₕ+1, Nₕ, Nₕ)
    tmp_Dstar   = zeros(lₕ+1, Nₕ, Nₕ)
    tmp_U       = zeros(lₕ+1, Nₕ, Nₕ)
    tmp_ϵ       = zeros(lₕ+1,Nₕ)
    tmp_ϵ_sort  = zeros((lₕ+1)*Nₕ)
    tmp_n       = zeros(lₕ+1, Nₕ)
    
    tmp_tn      = 0.0     

    CacheODA(A, M₀, M₋₁, M₋₂, Hfix, tmp_H, tmp_Dstar, tmp_D, tmp_U, tmp_ϵ, tmp_ϵ_sort, tmp_n, tmp_tn)
end


function performstep!(method::ODA, solver::KhonShamSolver)

    @unpack tmp_D, tmp_U, tmp_ϵ, tmp_n = solver.cache

    # STEP 1 : Resolution of the generalized eigenvalue problem to find atomic orbitals and corresonding energies
    find_orbital!(solver.discretization, solver)

    # STEP 2 : Build the n matrix using the Aufbau principle
    aufbau!(solver)

    # STEP 3 : Build the density related matrix
    update_density!(method, solver)

    # Registering into solver
    solver.D .= tmp_D
    solver.U .= tmp_U
    solver.ϵ .= tmp_ϵ
    solver.n .= tmp_n
end

stopping_criteria(m::ODA, solver::KhonShamSolver) = stopping_criteria(m, solver.D, solver.Dprev)
stopping_criteria(::ODA, D, Dprev) = norm(D .- Dprev)

function aufbau!(solver::KhonShamSolver)

    @unpack N = solver.model
    @unpack tmp_ϵ, tmp_n = solver.cache

    (d1,_) = size(tmp_ϵ)
    index_ϵ_sort = sortperm(tmp_ϵ, dims = 2)
    count = N
    i = firstindex(index_ϵ_sort)
    while count > 0
        l = div(index_ϵ_sort[i], d1)
        if count - (2*l + 1) ≥ 0
            for k in -l:l
                tmp_n[l+1,k+l+1] = 2
            end
            count = count - (2*l + 1)
        else
            for k in -l:l
                tmp_n[l+1,k+l+1] = 2/(2*l+1)
            end
            count = 0
            break
        end
        i+=1
    end
end


function find_orbital!(discretization::KohnShamSphericalDiscretization, solver::KhonShamSolver)

    @unpack lₕ = discretization
    @unpack M₀, M₋₁, M₋₂, Hfix, tmp_H, tmp_ϵ = solver.cache

    # STEP 1 : Find Hartree term 
    #build_hartree!(discretization, matrix) 

    # STEP 2 : compute an approximation of the exchange correlation term
    #build_exchange_corr!(discretization, matrix)

    for l ∈ 0:lₕ
        tmp_H .= Hfix[l+1]
        tmp_ϵ[l+1,:], tmp_U[l+1,:,:] = solve_generalized_eigenvalue_problem(tmp_H, M₀)
    end

end


function update_density(::ODA, solver::KhonShamSolver)

    # STEP 7 : update this matrix with a convex approach
    # some optimisation to find a good tₙ
    tmp_tₙ = 0.5
    @. tmp_D = tmp_tₙ * tmp_Dstar + (1 - tmp_tₙ) * solver.Dprev
end

