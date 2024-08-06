function aufbau!(solver::KhonShamSolver)
    @unpack N = solver.model
    @unpack tmp_ϵ, tmp_index_sort, tmp_n = solver.discretization.cache
    tmp_n .= zero(tmp_n)
    tmp_index_sort = aufbau!(tmp_n, tmp_ϵ, N; tol = solver.opts.degen_tol)
end

function aufbau!(n, ϵ, N; tol = eps(eltype(ϵ)))
    ϵ_copy = copy(ϵ)
    _l,_n = size(ϵ_copy)
    ϵ_vect = vec(ϵ_copy)
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
        if remain - 2 * degen ≥ 0
            for i in A
                n[i] = 2 * degen_matrix[i]
            end
            remain -= 2 * degen 
        else
            if length(A) == 1
                # First case, if no degeneracy
                n[first(A)] = remain
            else
                # Second case, if degeneracy
                @error "There is accidental degeneracy but no implementation for this case for the moment."
            end
            break
        end
    end
    index_sort
end