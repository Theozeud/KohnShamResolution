# OPTIMAL DAMPLING ALGORITHM

mutable struct ODA{T<:Real} <: RCAMethod 
    t::T
    iter::Int
    function ODA(t::Real, iter::Int = 1)
        new{typeof(t)}(t,iter)
    end
end

function update_density!(solver::KohnShamSolver, m::ODA)
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


function find_minima_oda(energy_kin0::Real, energy_kin1::Real, 
                         energy_cou0::Real, energy_cou1::Real, 
                         energy_har0::Real, energy_har1::Real, 
                         energy_har01::Real, energy_har10::Real,
                         D0::AbstractArray{<:Real}, D1::AbstractArray{<:Real},
                         model::AbstractDFTModel, discretization::KohnShamDiscretization)

    # DEFINE THE OBJECTIV FUNCTION
    function f(t::Float64)        
        (t * (energy_kin0 + energy_cou0) + (1-t) * (energy_kin1 + energy_cou1) + 
        t^2*energy_har0 + (1-t)^2 * energy_har1 + t*(1-t) * (energy_har01 + energy_har10) +
        compute_exchangecorrelation_energy(discretization, model, t*D0 .+ (1-t) * D1))
    end

    # PERFORM THE OPTIMISATION THROUGHT THE GOLDEN SECTION's METHOD
    res = optimize(f, 0.0, 1.0, GoldenSection())

    res.minimizer, res.minimum
end

function update_energy!(solver::KohnShamSolver, t::Real, Dt::AbstractArray{<:Real},
                        energy_kin0::Real, energy_kin1::Real, 
                        energy_cou0::Real, energy_cou1::Real, 
                        energy_har0::Real, energy_har1::Real, 
                        energy_har01::Real, energy_har10::Real)

    solver.energy_kin = t*energy_kin0 + (1-t)*energy_kin1
    solver.energy_cou = t*energy_cou0 + (1-t)*energy_cou1
    solver.energy_har = t^2*energy_har0 + (1-t)^2*energy_har1 + t*(1-t) * (energy_har01 + energy_har10)
    solver.energy_exc = compute_exchangecorrelation_energy(solver.discretization, solver.model, Dt)
    nothing
end


function update_energy!(solver::KohnShamSolver, t::Real, Dt::AbstractArray{<:Real},
                        energy_kin0::Real, energy_kin1::Real, 
                        energy_cou0::Real, energy_cou1::Real)

    solver.energy_kin = t*energy_kin0 + (1-t)*energy_kin1
    solver.energy_cou = t*energy_cou0 + (1-t)*energy_cou1
    solver.energy_har = compute_hartree_energy(solver.discretization, Dt)
    solver.energy_exc = compute_exchangecorrelation_energy(solver.discretization, solver.model, Dt)
    nothing
end