function aufbau!(solver::KhonShamSolver)
    aufbau!(solver.discretization, solver.model.N; tol = solver.opts.degen_tol)
end

function aufbau!(discretization::LDADiscretization, N; tol = eps(eltype(ϵ)))
    @unpack tmp_ϵ, tmp_index_sort, tmp_n = discretization.tmp_cache
    tmp_index_sort .= sortperm(vec(tmp_ϵ))
    function degeneracy(idx)
        l = rem(idx - 1, size(tmp_ϵ, 1)) 
        return 4 * l + 2
    end
    remain = N
    idx = 1
    while remain > 0 && idx ≤ length(tmp_index_sort)
        indices_degen = [tmp_index_sort[idx]]  
        idx += 1
        while idx ≤ length(tmp_index_sort) && abs(tmp_ϵ[tmp_index_sort[idx]] - tmp_ϵ[first(indices_degen)]) < tol
            push!(indices_degen, tmp_index_sort[idx])
            idx += 1
        end
        # Count total degeneracy
        total_degen = sum(degeneracy(i) for i in indices_degen)

        # Electron distribution
        if remain - total_degen≥ 0
            for i in indices_degen
                tmp_n[i] = degeneracy(i)
            end
            remain -= total_degen
        else
            if length(indices_degen) == 1
                # First case, if no degeneracy
                tmp_n[first(indices_degen)] = remain
            else
                # Second case, if degeneracy
                @error "There is accidental degeneracy but no implementation for this case for the moment."
            end
            remain = zero(remain)
        end
    end
    nothing
end


function aufbau!(discretization::LSDADiscretization, N; tol = eps(eltype(ϵ)))
    @unpack tmp_ϵ, tmp_index_sort, tmp_n = discretization.tmp_cache
    tmp_index_sort .= sortperm(vec(tmp_ϵ))
    function degeneracy(idx)
        l = rem(idx - 1, size(tmp_ϵ, 1)) 
        return 2 * l + 1
    end
    remain = N
    idx = 1
    while remain > 0 && idx ≤ length(tmp_index_sort)
        indices_degen = [tmp_index_sort[idx]]  
        idx += 1
        while idx ≤ length(tmp_index_sort) && abs(tmp_ϵ[tmp_index_sort[idx]] - tmp_ϵ[first(indices_degen)]) < tol
            push!(indices_degen, tmp_index_sort[idx])
            idx += 1
        end
        # Count total degeneracy
        total_degen = sum(degeneracy(i) for i in indices_degen)

        # Electron distribution
        if remain - total_degen≥ 0
            for i in indices_degen
                tmp_n[i] = degeneracy(i)
            end
            remain -= total_degen
        else
            tmp_n[1] = remain
            #=
            for i in indices_degen
                tmp_n[i] = remain/length(indices_degen)
            end
            =#
            #=
            if length(indices_degen) == 1
                # First case, if no degeneracy
                tmp_n[first(indices_degen)] = remain
            else
                # Second case, if degeneracy
                @error "There is accidental degeneracy but no implementation for this case for the moment."
            end
            =#
            remain = zero(remain)
        end
    end
    nothing
end