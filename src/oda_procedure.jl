function find_minima_oda(energy_kin0::Real, energy_kin1::Real, 
                         energy_cou0::Real, energy_cou1::Real, 
                         energy_har0::Real, energy_har1::Real, energy_har01::Real,
                         D0::AbstractArray{<:Real}, D1::AbstractArray{<:Real},
                         discretization::KohnShamDiscretization, model::KohnShamExtended)
    println(energy_kin0 + energy_cou0)
    println(energy_kin1 + energy_cou1)
    println(energy_har0)
    println(energy_har1)
    println(energy_har01)
    # DEFINE THE OBJECTIV FUNCTION
    function f(t::Float64)        
        (t * (energy_kin0 + energy_cou0) + (1-t) * (energy_kin1 + energy_cou1) + 
        t^2*energy_har0 + (1-t)^2 * energy_har1 + t*(1-t) * energy_har01 +
        compute_exchangecorrelation_energy(discretization, model, t*D0 .+ (1-t) * D1))
    end

    # PERFORM THE OPTIMISATION THROUGHT THE GOLDEN SECTION's METHOD
    res = optimize(f, 0.0, 1.0, GoldenSection())
    res2 = optimize(f, res.minimizer, 1.0, GoldenSection())
    res2.minimizer
end