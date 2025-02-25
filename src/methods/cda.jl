# CONSTANT DAMPLING ALGORITHM

struct CDA{T} <: RCAMethod
    t::T
    function CDA(t::Real)
        @assert 0 ≤ t ≤ 1
        new{typeof(t)}(t)
    end
end

function update_density!(m::CDA, solver::KhonShamSolver)
    @unpack D, Dprev, energy_kin, energy_cou,
            energy_kin_prev, energy_cou_prev  = solver

    if solver. niter > 0
        # UPDATE THE DENSITY
        @. D = m.t * D + (1 - m.t) * Dprev

        # UPDATE THE ENERGIES
        update_energy!( solver, m.t, D,
                        energy_kin, energy_kin_prev, 
                        energy_cou, energy_cou_prev)

        solver.energy_har = compute_hartree_energy(solver.discretization, D)
        solver.energy = m.t * solver.energy + (1-m.t) * solver.energy_prev
    end
    nothing
end