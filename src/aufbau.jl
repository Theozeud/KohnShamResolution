
# AUFBAU PRINCIPLE
function aufbau!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack ϵ, n, D, tmpD, U, model, opts = solver
    @unpack N = model
    @unpack degen_tol = opts
    @unpack tmp_index_sort = discretization.tmp_cache

    n  .= zero(n)

    # RETURN THE MAXIMAL NUMBER OF ELECTRONS THAT CAN BE PUT AT THE LAYER CORRESPONDING TO THE INDEX idx
    function degeneracy(idx)
        l = rem(idx - 1, size(ϵ, 1)) 
        return 4 * l + 2
    end

    # COLLECT THE INDICES OF THE SORTED ORBITAL ENERGIES
    tmp_index_sort .= sortperm(vec(ϵ))
   
    # LOOP TO FILL THE ORBITALS
    remain = N
    idx = 1
    while remain > 0 && idx ≤ length(tmp_index_sort)

        # FIND ALL THE ORBITALS WITH THE SAME ENERGY
        indices_degen = [tmp_index_sort[idx]]  
        idx += 1
        while idx ≤ length(tmp_index_sort) && abs(ϵ[tmp_index_sort[idx]] - ϵ[first(indices_degen)]) < degen_tol
            push!(indices_degen, tmp_index_sort[idx])
            idx += 1
        end

        # COUNT TOTAL DEGENERACY
        total_degen = sum(degeneracy(i) for i in indices_degen)

        # ELECTRON DISTRIBUTIONS
        if remain - total_degen ≥ 0
            # IN THIS CASE WE CAN FILL ALL THE LAYERS WITH THE SAME ORBITAL ENERGY
            for i in indices_degen
                n[i] = degeneracy(i)
            end
            remain -= total_degen
        else
            println(indices_degen)
            # IN THIS CASE WE NEED TO FILL THE LAYERS WITH THE OPTIMAL REPARTITION
            if length(indices_degen) == 1
                # NO DEGENERACY
                println("Degen 1")
                n[first(indices_degen)] = remain

            elseif length(indices_degen) == 2
                # DEGENERACY OF ORDER 2
                println("Degen 2")
                # COMPUTE ONE EXTREMA 
                degen_first = degeneracy(indices_degen[1])
                if degen_first ≥ remain
                    n[indices_degen[1]] = remain
                    n[indices_degen[2]] = zero(remain)
                    n1_0 = remain
                    n2_0 = zero(remain)
                else
                    n[indices_degen[1]] = degen_first
                    n[indices_degen[2]] = remain - degen_first
                    n1_0 = degen_first
                    n2_0 = remain - degen_first
                end
                
                density_matrix!(discretization, U, n, D)

                energy_kin0 = compute_kinetic_energy(discretization, U, n)
                energy_cou0 = compute_coulomb_energy(discretization, U, n)
                energy_har0 = compute_hartree_energy(discretization, D)

                # COMPUTE THE OTHER EXTREMA
                degen_second = degeneracy(indices_degen[2])
                if degen_second ≥ remain
                    n[indices_degen[1]] = zero(remain)
                    n[indices_degen[2]] = remain
                    n1_1 = zero(remain)
                    n2_1 = remain
                else
                    n[indices_degen[1]] = degen_first
                    n[indices_degen[2]] = remain - degen_first
                    n1_1 = remain - degen_second
                    n2_1 = degen_second
                end
                
                density_matrix!(discretization, U, n, tmpD)

                energy_kin1 = compute_kinetic_energy(discretization, U, n)
                energy_cou1 = compute_coulomb_energy(discretization, U, n)
                energy_har1 = compute_hartree_energy(discretization, tmpD)

                energy_har01 = compute_hartree_mix_energy(discretization, D, tmpD)

                # FIND THE OPTIMUM OCCUPATION
                t = find_minima_oda(energy_kin0, energy_kin1, 
                                    energy_cou0, energy_cou1, 
                                    energy_har0, energy_har1, energy_har01,
                                    D, tmpD, discretization, model)
                println(t)
                n[indices_degen[1]] = t * n1_0 + (1-t) * n1_1
                n[indices_degen[2]] = t * n2_0 + (1-t) * n2_1
                
            else
                @error("This case of degeneracy is not coded.")
            end
            remain = zero(remain)
        end
    end
    nothing
end


function aufbau!(discretization::LSDADiscretization, ϵ, n, N; tol = eps(eltype(ϵ)))
    @unpack tmp_index_sort = discretization.tmp_cache
    tmp_index_sort .= sortperm(vec(ϵ))
    function degeneracy(idx)
        l = rem(idx - 1, size(ϵ, 1)) 
        return 2 * l + 1
    end
    remain = N
    idx = 1
    while remain > 0 && idx ≤ length(tmp_index_sort)
        indices_degen = [tmp_index_sort[idx]]  
        idx += 1
        while idx ≤ length(tmp_index_sort) && abs(ϵ[tmp_index_sort[idx]] - ϵ[first(indices_degen)]) < tol
            push!(indices_degen, tmp_index_sort[idx])
            idx += 1
        end
        # Count total degeneracy
        total_degen = sum(degeneracy(i) for i in indices_degen)

        # Electron distribution
        if remain - total_degen≥ 0
            for i in indices_degen
                n[i] = degeneracy(i)
            end
            remain -= total_degen
        else
            n[1] = remain
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