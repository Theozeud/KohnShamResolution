# Constant Dampling Algotithm
struct CDA{T} <: SCFMethod
    t::T
    function CDA(t::Real)
        @assert 0 ≤ t ≤ 1
        new{typeof(t)}(t)
    end
end

function update_density!(m::CDA, solver::KhonShamSolver)
    @unpack D, Dprev = solver
    @unpack t = m
    @. D = t * D + (1 - t) * Dprev
    nothing
end

# OPTIMAL DAMPLING ALGORITHM

mutable struct ODA{T<:Real} <: SCFMethod 
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

    if solver.niter < m.iter
        D .= m.t * D + (1 - m.t) * Dprev
        return nothing
    end
    
    energy_har01 = compute_hartree_mix_energy(discretization, D, Dprev)

    # FIND THE OPTIMUM OCCUPATION
    m.t = find_minima_oda(energy_kin, energy_kin_prev, 
                          energy_cou, energy_cou_prev, 
                          energy_har, energy_har_prev, energy_har01,
                          D, Dprev, discretization, model)
    @show m.t

    # UPDATE THE DENSITY
    D .= m.t * D + (1 - m.t) * Dprev

    # UPDATE THE ENERGY
    solver.energy_kin = m.t * energy_kin + (1-m.t) * energy_kin_prev
    solver.energy_cou = m.t * energy_cou + (1-m.t) * energy_cou_prev
    solver.energy_har = elT(0.5)*(m.t^2 * energy_har + (1-m.t)^2*energy_har_prev) 
                        + m.t*(1-m.t)*energy_har01
    solver.energy_exc = compute_exchangecorrelation_energy(discretization, model, D)
    nothing
end
