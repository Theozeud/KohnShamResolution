# OPTIMAL DAMPLING ALGORITHM

mutable struct ODA{T<:Real} <: RCAMethod 
    t::T
    iter::Int
    function ODA(t::Real, iter::Int = 1)
        new{typeof(t)}(t,iter)
    end
end


function update_density!(m::ODA, solver::KhonShamSolver)
    @unpack D, Dprev, discretization, model,
            energy_kin, energy_cou, energy_har,
            energy_kin_prev, energy_cou_prev, energy_har_prev = solver
    @unpack elT = discretization

    if solver.niter > 0
        
        if solver.niter < m.iter
            D .= m.t * D + (1 - m.t) * Dprev

            # UPDATE THE ENERGIES
            update_energy!( solver, m.t, D,
                            energy_kin, energy_kin_prev, 
                            energy_cou, energy_cou_prev)

            solver.energy_har = compute_hartree_energy(discretization, D)
            solver.energy = m.t * solver.energy + (1-m.t) * solver.energy_prev
            return nothing
        end
        
        energy_har01 = compute_hartree_mix_energy(discretization, D, Dprev)
        energy_har10 = compute_hartree_mix_energy(discretization, Dprev, D)

        # FIND THE OPTIMUM OCCUPATION
        m.t, energy = find_minima_oda(energy_kin, energy_kin_prev, 
                                      energy_cou, energy_cou_prev, 
                                      energy_har, energy_har_prev, 
                                      energy_har01, energy_har10,
                                      D, Dprev, model, discretization)
        @show m.t

        # UPDATE THE DENSITY
        D .= m.t * D + (1 - m.t) * Dprev

        # UPDATE THE ENERGIES
        update_energy!( solver, m.t, D,
                        energy_kin, energy_kin_prev, 
                        energy_cou, energy_cou_prev, 
                        energy_har, energy_har_prev, 
                        energy_har01, energy_har10)

        solver.energy = energy
    end
    nothing
end