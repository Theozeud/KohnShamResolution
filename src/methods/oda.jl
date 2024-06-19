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
    tmp_index_sort
    tmp_n    
    tmp_tn
end

function init_cache(::ODA, model::AbstractDFTModel, discretization::KohnShamDiscretization)

    @unpack lₕ, Nₕ, basis, mesh, Rmin, Rmax = discretization

    # Set the type of number as the one of the discretization basis
    T = bottom_type(discretization.basis)

    # Init base matrices
    @assert length(basis) == Nₕ
    
    deriv_basis = deriv(basis)

    A   = mass_matrix(deriv_basis, Rmin, Rmax)
    M₀  = mass_matrix(basis, Rmin, Rmax)
    M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
    M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

    # Creation of the fix part of the hamiltonian   
    Kin =  zeros(T, lₕ+1, Nₕ, Nₕ)
    Coulomb =  zeros(T, lₕ+1, Nₕ, Nₕ)
    build_kinetic!(discretization, Kin, A, M₋₂)
    build_coulomb!(discretization, Coulomb, model, M₋₁)
    Hfix = Kin + Coulomb

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, Nₕ, Nₕ)
    tmp_D           = zero_piecewiselaurantpolynomial(mesh, T)
    tmp_Dstar       = zero_piecewiselaurantpolynomial(mesh, T)
    tmp_U           = zeros(T, lₕ+1, Nₕ, Nₕ)
    tmp_exc         = zeros(T, Nₕ, Nₕ)
    tmp_Hartree     = zeros(T, Nₕ, Nₕ)
    tmp_ϵ           = zeros(T, lₕ+1, Nₕ)
    tmp_index_sort  = zeros(Int, Nₕ*(lₕ+1))
    tmp_n           = zeros(T, lₕ+1, Nₕ)
        
    tmp_tn      = T(0)     

    CacheODA(A, M₀, M₋₁, M₋₂, Hfix, tmp_H, tmp_Dstar, tmp_D, tmp_U, tmp_Hartree, tmp_exc, tmp_ϵ, tmp_index_sort, tmp_n, tmp_tn)
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
stopping_criteria(::ODA, D, Dprev, points) = sqrt((points[end] - points[begin])*sum([abs(D(x)- Dprev(x))^2 for x ∈ points[begin:end-1]])/(length(points)-1))

function find_orbital!(discretization::KohnShamSphericalDiscretization, solver::KhonShamSolver)

    @unpack lₕ = discretization
    @unpack M₀, Hfix, tmp_H, tmp_U, tmp_Hartree, tmp_exc, tmp_ϵ = solver.cache
    @unpack exc = solver.model
    @unpack Dprev = solver
    @unpack quad_method, quad_reltol, quad_abstol, hartree = solver.opts

    # STEP 1 : Compute Hartree term 
    if hartree
        build_hartree!(discretization, tmp_Hartree, Dprev)
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
    @unpack tmp_ϵ, tmp_index_sort, tmp_n = solver.cache
    @unpack lₕ, Nₕ = solver.discretization
    tmp_index_sort = aufbau!(tmp_n, tmp_ϵ, N; tol = solver.opts.degen_tol)
end

function aufbau!(n, ϵ, N; tol = eps(eltype(ϵ)))
    _l,_n = size(ϵ)
    ϵ_vect = vec(ϵ)
    index_sort = sortperm(ϵ_vect)
    degen_matrix = reduce(hcat, [[2*l + 1 for l ∈ 0:_l-1] for i ∈ 1:_n])
    remain = N
    idx = 1
    while remain > 0 && idx < length(ϵ) + 1
        A = Int[]  #Stock all index corresponding to a degenerancy
        ϵ_cur = ϵ_vect[index_sort[idx]]
        push!(A, index_sort[idx])
        idx += 1
        while abs(ϵ[index_sort[idx]] - ϵ_cur) < tol
            push!(A, index_sort[idx])
            idx += 1
        end
        # Count total degeneracy
        degen = sum(degen_matrix[i] for i in A)
        # See what to do depending on the case
        if remain - degen ≥ 0
            for i in A
                n[i] = 2 * degen_matrix[i]
            end
            remain -= degen
        else
            if length(A) == 1
                # First case, if no degeneracy
                n[first(A)] = 2 * remain
            else
                # Second case, if degeneracy
                @error "There is accidental degeneracy but no implementation for this case for the moment."
            end
            break
        end
    end
    index_sort
end

function update_density!(::ODA, solver::KhonShamSolver)

    @unpack Dprev = solver
    @unpack tmp_D, tmp_Dstar, tmp_U, tmp_n, tmp_tn = solver.cache

    tmp_Dstar = build_density!(solver.discretization, tmp_Dstar, tmp_U, tmp_n)
    tmp_tn = 1.0
    tmp_D = tmp_tn * tmp_Dstar + (1 - tmp_tn) * Dprev
end