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
    function ODA(t::Real)
        new{typeof(t)}(t)
    end
end


function update_density!(m::ODA, solver::KhonShamSolver)
    @unpack D, Dprev = solver
    @unpack energy_kin, energy_cou, energy_har,
    energy_kin_prev, energy_cou_prev, energy_har_prev = solver
    if iszero(solver.niter)
        D .= m.t * D + (1 - m.t) * Dprev
        return nothing
    end
    if !isthereExchangeCorrelation(solver.model)
        energy_har_mix = compute_hartree_mix_energy(solver.discretization, solver)
        t2 = energy_har - 2*energy_har_mix + energy_har_prev 
        t1 = energy_kin - energy_kin_prev - 2 * energy_har_prev
        t0 = energy_kin_prev + energy_har_prev
        if t2 > 0
            _argmin = -t1/(2*t2)
            m.t = min(max(0.0,_argmin),1.0)
        elseif t2 < 0
            if t2*0.1^2 + t1*0.1 > t2 + t1
                m.t = 1.0
            else
                m.t = 0.1
            end
        else
            if t1 > 0
                m.t = 0.1
            elseif t1 < 0
                m.t = 1.0
            else
                m.t = 1.0
            end
        end
    else

    end
    @show m.t
    D .= m.t * D + (1 - m.t) * Dprev
    nothing
end
