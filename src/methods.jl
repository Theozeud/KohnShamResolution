# Constant Dampling Algotithm
struct CDA{T} <: SCFMethod
    t::T
    function CDA(t::Real)
        @assert 0 ≤ t ≤ 1
        new{typeof(t)}(t)
    end
end

function update_density!(m::CDA, solver::KhonShamSolver)
    @unpack tmp_D, tmp_Dstar = solver.discretization.tmp_cache
    tmp_D .= m.t * tmp_Dstar + (1 - m.t) * solver.Dprev
    nothing
end

# Optimal Dampling Algotithm


