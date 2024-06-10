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

    @unpack lₕ, Nₕ, basis, mesh = discretization

    # Init base matrices
    @assert length(basis) == Nₕ
    
    deriv_basis = deriv(basis)

    A   = mass_matrix(deriv_basis, mesh[begin], mesh[end])
    M₀  = mass_matrix(basis, mesh[begin], mesh[end])
    M₋₁ = weight_mass_matrix(basis, -1, mesh[begin], mesh[end])
    M₋₂ = weight_mass_matrix(basis, -2, mesh[begin], mesh[end])

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
    tmp_ϵ       = zeros(lₕ+1, Nₕ)
    tmp_ϵ_sort  = Int[]
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
    @unpack tmp_ϵ, tmp_ϵ_sort, tmp_n = solver.cache
    @unpack lₕ, Nₕ = solver.discretization

    # splat the vector of eigenvalue
    vector_ϵ = [tmp_ϵ...]

    # Create the degeneracy matrix
    degen_matrix = reduce(hcat, [[2*l + 1 for l ∈ 0:lₕ] for i ∈ 1:Nₕ])

    remain = N
    idx = 1
    while remain > 0 && idx < (lₕ+1)*Nₕ
        
        # Find the next lowest eigenvalue(s) (may be several in case of degeneracy)
        a = first(vector_ϵ)
        A = [1]
        for j ∈ eachindex(vector_ϵ)
            if vector_ϵ[j] < a 
                a = vector_ϵ[j]
                A = [j]
            elseif vector_ϵ[j] == a 
                push!(A,j)
            end
        end

        # Remove this eigenvalue
        for i in A
            deleteat!(vector_ϵ, i)
        end
        
        # Count total degeneracy
        degen = sum(degen_matrix[i] for i in A)
        
        # See what to do depending on the case
        if remain - degen ≥ 0
            for i in A
                tmp_n[i] = 2 * degen[i]
            end
            remain -= degen
        else
            if length(A) == 1
                # First case, if no degeneracy
                tmp_n[first(A)] = 2 * remain
            else
                # Second case, if degeneracy
                @error "There is degeneracy but no implementation for this case for the moment."
            end
            break
        end
        idx += length(A)
    end
end


function find_orbital!(discretization::KohnShamSphericalDiscretization, solver::KhonShamSolver)

    @unpack lₕ = discretization
    @unpack M₀, M₋₁, M₋₂, Hfix, tmp_H, tmp_U, tmp_ϵ = solver.cache

    # STEP 1 : Find Hartree term 
    #build_hartree!(discretization, matrix) 

    # STEP 2 : compute an approximation of the exchange correlation term
    #build_exchange_corr!(discretization, matrix)

    # STEP 3 : Solve the generalized eigenvalue problem for each section l
    for l ∈ 0:lₕ
        # building the hamiltonian of the lᵗʰ section
        tmp_H .= Hfix[l+1]
        # solving
        tmp_ϵ[l+1,:], tmp_U[l+1,:,:] = solve_generalized_eigenvalue_problem(tmp_H, M₀)
    end
end


function update_density!(::ODA, solver::KhonShamSolver)

    @unpack Dprev = solver
    @unpack tmp_D, tmp_Dstar, tmp_tn = solver.cache

    tmp_tn = 0.5
    @. tmp_D = tmp_tn * tmp_Dstar + (1 - tmp_tn) * solver.Dprev
end

