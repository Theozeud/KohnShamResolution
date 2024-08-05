abstract type ODA <: AbstractKohnShamResolutionMethod end

struct ConstantODA{T} <: ODA
    t::T
    function ConstantODA(t::Real = 1)
        @assert 0≤t≤ 1
        new{typeof(t)}(t)
    end
end

#= ODA not coded
function update_density!(::ODA, solver::KhonShamSolver)

    @unpack Dprev = solver
    @unpack tmp_D, tmp_Dstar, tmp_U, tmp_n, tmp_tn = solver.discretization.cache

    tmp_Dstar = build_density!(solver.discretization, tmp_Dstar, tmp_U, tmp_n)
    tmp_tn = 0.5
    tmp_D = tmp_tn * tmp_Dstar + (1 - tmp_tn) * Dprev
end
=#

function update_density!(m::ConstantODA, solver::KhonShamSolver)
    @unpack tmp_D, tmp_Dstar = solver.discretization.cache
    tmp_D .= m.t * tmp_Dstar + (1 - m.t) * solver.Dprev
    nothing
end