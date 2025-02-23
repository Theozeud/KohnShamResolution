# RETURN THE MAXIMAL NUMBER OF ELECTRONS THAT CAN BE PUT AT THE LAYER CORRESPONDING TO THE INDEX idx

function convert_index(discretization::LDADiscretization, idx::Int)
    @unpack lₕ  = discretization
    l = rem(idx - 1, lₕ+1)
    k = div(idx-1,lₕ+1)+1
    return (l,k)
end

function degeneracy(discretization::LDADiscretization, idx::Int)
    l,_ = convert_index(discretization::LDADiscretization, idx)
    return 4 * l + 2
end

function degeneracy(::LSDADiscretization, idx::Int)
    l = rem(idx - 1, size(ϵ, 1)) 
    return 2 * l + 1
end

# AUFBAU PRINCIPLE
function aufbau!(discretization::LDADiscretization, solver::KhonShamSolver)

    @unpack ϵ, n, D, tmpD, U, model, opts = solver
    @unpack tmp_index_sort = discretization.tmp_cache

    # INIT OCCUPATION NUMBER 
    n  .= zero(n)

    # COLLECT THE INDICES OF THE SORTED ORBITAL ENERGIES
    tmp_index_sort .= sortperm(vec(ϵ))
   
    # LOOP TO FILL THE ORBITALS
    remain = model.N
    idx = 1
    while remain > 0 && idx ≤ length(tmp_index_sort)

        # FIND ALL THE ORBITALS WITH THE SAME ENERGY
        indices_degen = [tmp_index_sort[idx]]  
        idx += 1
        while idx ≤ length(tmp_index_sort) && abs(ϵ[tmp_index_sort[idx]] - ϵ[first(indices_degen)]) < opts.degen_tol && length(indices_degen) <2
            push!(indices_degen, tmp_index_sort[idx])
            idx += 1
        end

        # COMPUTE DEGENERACY
        degen = zeros(Int,length(indices_degen))
        for i ∈ eachindex(indices_degen)
            degen[i] = degeneracy(discretization, indices_degen[i])
        end
        total_degen = sum(degen)

        # ELECTRON DISTRIBUTIONS
        if remain - total_degen ≥ 0
            # IN THIS CASE WE CAN FILL ALL THE LAYERS WITH THE SAME ORBITAL ENERGY

            for i in eachindex(indices_degen)
                n[indices_degen[i]] = degen[i]
                normalization!(discretization, solver, convert_index(discretization,i)...) 
            end
            remain -= total_degen

        elseif length(indices_degen) == 1
            # IN THIS CASE, WE FILL THE ORBITAL WITH THE REMAIN ELECTRONS

            n[first(indices_degen)] = remain
            normalization!(discretization, solver, convert_index(discretization,first(indices_degen))...)
            break

        elseif length(indices_degen) == 2
            # IN THIS CASE WE NEED TO FILL THE LAYERS WITH THE OPTIMAL REPARTITION
            println("Degen 2")
            solver.flag_degen = true

            # NORMALIZATION OF COEFFICIENTS OF ORBITAL
            for i ∈ indices_degen
                normalization!(discretization, solver, convert_index(discretization,i)...)
            end

            # COMPUTE ENERGIES FOR ONE EXTREMA 
            degen_first = degen[1]
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

            # COMPUTE ENERGIES FOR THE OTHER EXTREMA
            degen_second = degen[2]
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
            energy_har10 = compute_hartree_mix_energy(discretization, tmpD, D)

            # FIND THE OPTIMUM OCCUPATION
            t, energy = find_minima_oda(energy_kin0, energy_kin1, 
                                        energy_cou0, energy_cou1, 
                                        energy_har0, energy_har1, 
                                        energy_har01, energy_har10,
                                        D, tmpD, model, discretization)
            println(t)

            # UPDATE THE OCCUPATION NUMBERS
            n[indices_degen[1]] = t * n1_0 + (1-t) * n1_1
            n[indices_degen[2]] = t * n2_0 + (1-t) * n2_1

            # UPDATE THE DENSITY
            @. D = t * D + (1 - t) * tmpD

            # UPDATE THE ENERGIES
            update_energy!( solver, t, D,
                            energy_kin0, energy_kin1, 
                            energy_cou0, energy_cou1, 
                            energy_har0, energy_har1, 
                            energy_har01, energy_har10)

            solver.energy = energy
            break
        else
            @error("This case of degeneracy is not coded.")
            break
        end

    end
    nothing
end


