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
    tmp_Hartree
    tmp_exc 
    tmp_ϵ
    tmp_ϵ_sort
    tmp_n    
    tmp_tn
end

function init_cache(::ODA, model::AbstractDFTModel, discretization::KohnShamDiscretization)

    @unpack lₕ, Nₕ, basis, mesh, Rmin, Rmax = discretization

    # Init base matrices
    @assert length(basis) == Nₕ
    
    deriv_basis = deriv(basis)

    A   = mass_matrix(deriv_basis, Rmin, Rmax)
    M₀  = mass_matrix(basis, Rmin, Rmax)
    M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
    M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

    # Creation of the fix part of the hamiltonian   
    Kin =  zeros(lₕ+1, Nₕ, Nₕ)
    Coulomb =  zeros(lₕ+1, Nₕ, Nₕ)
    build_kinetic!(discretization, Kin, A, M₋₂)
    build_coulomb!(discretization, Coulomb, model, M₋₁)
    Hfix = Kin + Coulomb

    # Initialization of array for temporary stockage of computations
    tmp_H       = zeros(Nₕ, Nₕ)

    tmp_D       = zero_piecewiselaurantpolynomial(mesh)
    tmp_Dstar   = zero_piecewiselaurantpolynomial(mesh)
    tmp_U       = zeros(lₕ+1, Nₕ, Nₕ)
    tmp_exc     = zeros(Nₕ, Nₕ)
    tmp_Hartree = zeros(Nₕ, Nₕ)
    tmp_ϵ       = zeros(lₕ+1, Nₕ)
    tmp_ϵ_sort  = Int[]
    tmp_n       = zeros(lₕ+1, Nₕ)
        
    tmp_tn      = 0.0     

    CacheODA(A, M₀, M₋₁, M₋₂, Hfix, tmp_H, tmp_Dstar, tmp_D, tmp_U, tmp_Hartree, tmp_exc, tmp_ϵ, tmp_ϵ_sort, tmp_n, tmp_tn)
end

function performstep!(method::ODA, solver::KhonShamSolver)

    # STEP 1 : Resolution of the generalized eigenvalue problem to find atomic orbitals and corresonding energies
    find_orbital!(solver.discretization, solver)

    # STEP 2 : Build the n matrix using the Aufbau principle
    aufbau!(solver)

    # STEP 3 : Build the density related matrix
    @unpack tmp_D, tmp_Dstar, tmp_U, tmp_Hartree, tmp_exc, tmp_ϵ, tmp_n = solver.cache

    tmp_D = update_density!(method, solver)

    # Registering into solver
    solver.D  = tmp_D
    solver.U .= tmp_U
    solver.ϵ .= tmp_ϵ
    solver.n .= tmp_n

    tmp_D        = zero(tmp_D)
    tmp_Dstar    = zero(tmp_Dstar)
    tmp_U       .= zero(tmp_U)
    tmp_exc     .= zero(tmp_exc)
    tmp_Hartree .= zero(tmp_Hartree)
    tmp_ϵ       .= zero(tmp_ϵ)
    tmp_n       .= zero(tmp_n)
end

stopping_criteria(m::ODA, solver::KhonShamSolver) = stopping_criteria(m, solver.D, solver.Dprev, solver.discretization.mesh)
stopping_criteria(::ODA, D, Dprev, points) = sum([D(x)- Dprev(x) for x ∈ points])/length(points)

function find_orbital!(discretization::KohnShamSphericalDiscretization, solver::KhonShamSolver)

    @unpack lₕ = discretization
    @unpack M₀, Hfix, tmp_H, tmp_U, tmp_Hartree, tmp_exc, tmp_ϵ = solver.cache
    @unpack exc = solver.model
    @unpack Dprev = solver
    @unpack quad_method, quad_reltol, quad_abstol, hartree = solver.opts

    # STEP 1 : Compute Hartree term 
    if hartree
        build_hartree(discretization, tmp_Hartree, Dprev)
    end

    # STEP 2 : Compute Exchange Correlation term
    build_exchange_corr!(discretization, tmp_exc, Dprev, exc; quad_method = quad_method, quad_reltol = quad_reltol, quad_abstol = quad_abstol)

    # STEP 3 : Solve the generalized eigenvalue problem for each section l
    for l ∈ 0:lₕ
        # building the hamiltonian of the lᵗʰ section
        @. tmp_H = Hfix[l+1,:,:] + tmp_exc + tmp_Hartree
        # solving
        tmp_ϵ[l+1,:], tmp_U[l+1,:,:] = solve_generalized_eigenvalue_problem(tmp_H, M₀)
    end
end

function aufbau!(solver::KhonShamSolver)

    @unpack N = solver.model
    @unpack tmp_ϵ, tmp_ϵ_sort, tmp_n = solver.cache
    @unpack lₕ, Nₕ = solver.discretization

    # sort tmp_ϵ to order in l
    sort!(tmp_ϵ; dims = 1)

    # splat the vector of eigenvalue
    vector_ϵ = [tmp_ϵ...]

    # Create the degeneracy matrix
    degen_matrix = reduce(hcat, [[2*l + 1 for l ∈ 0:lₕ] for i ∈ 1:Nₕ])

    remain = N
    idx = 1
    B = Int[]
    C = collect(1:(lₕ+1)*Nₕ)
    while remain > 0 && idx < (lₕ+1)*Nₕ
        
        # Find the next lowest eigenvalue(s) (may be several in case of degeneracy)
        A = Int[]
        a = isempty(B) ? first(vector_ϵ) : vector_ϵ[first(C)]
        for j ∈ eachindex(vector_ϵ)
            if j ∉ B
                if vector_ϵ[j] < a 
                    a = vector_ϵ[j]
                    A = [j]
                elseif vector_ϵ[j] == a 
                    push!(A,j)
                end
            end
        end
        
        # Transfer to B the index of
        push!(B, A...)
        for i in A
            remove!(C, i)
        end

        # Count total degeneracy
        degen = sum(degen_matrix[i] for i in A)
        
        # See what to do depending on the case
        if remain - degen ≥ 0
            for i in A
                tmp_n[i] = 2 * degen_matrix[i]
            end
            remain -= degen
        else
            if length(A) == 1
                # First case, if no degeneracy
                tmp_n[first(A)] = 2 * remain
            else
                # Second case, if degeneracy
                @error "There is accidental degeneracy but no implementation for this case for the moment."
            end
            break
        end
        idx += length(A)
    end
end

function update_density!(::ODA, solver::KhonShamSolver)

    @unpack Dprev = solver
    @unpack tmp_D, tmp_Dstar, tmp_U, tmp_n, tmp_tn = solver.cache

    tmp_Dstar = build_density!(solver.discretization, tmp_Dstar, tmp_U, tmp_n)
    tmp_tn = 0.8
    tmp_D = tmp_tn * tmp_Dstar + (1 - tmp_tn) * Dprev
end