function find_orbital!(discretization::KohnShamRadialDiscretization, solver::KhonShamSolver)

    @unpack lₕ = discretization
    @unpack M₀, Hfix, Hartree, Exc, tmp_H, tmp_U, tmp_ϵ = discretization.cache
    @unpack Dprev = solver
    @unpack hartree, potential = solver.opts

    # STEP 1 : Compute Hartree term 
    if !iszero(hartree)
        hartree_matrix_pde!(discretization, Dprev)
        @. Hartree = hartree * Hartree
    end

    # STEP 2 : Compute Exchange Correlation term
    if isthereExchangeCorrelation(solver.model)
        exchange_corr_matrix!(discretization, solver.model)
    end

    # STEP 3 : Solve the generalized eigenvalue problem for each section l
    for l ∈ 0:lₕ
        # building the hamiltonian of the lᵗʰ section
        @. tmp_H = Hfix[l+1,:,:] + Exc + Hartree
        # solving
        tmp_ϵ[l+1,:], tmp_U[l+1,:,:] = solve_generalized_eigenvalue_problem(tmp_H, M₀)
    end
end