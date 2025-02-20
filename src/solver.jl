struct SolverOptions{T}
    scftol::T                       # SCF tolerance
    maxiter::Int                    # Maximum of iteration done
    quad_method                     # Method to use for quadrature of integrals
    quad_reltol::T                  # Relativ tolerance for the quadrature of integrals
    quad_abstol::T                  # Absolute tolerance for the quadrature of integrals
    hartree::T                      # Coefficient multiply to the Hartree Matrix : 
                                    # 0 -> no hartree term, 1-> full hartree term
    degen_tol::T    
end


mutable struct KhonShamSolver{  discretizationType <: KohnShamDiscretization,
                                modelType <: AbstractDFTModel,
                                methodType <: SCFMethod,
                                optsType <: SolverOptions,
                                logbookType <: LogBook,
                                dataType <: Real,
                                densityType <: AbstractArray,
                                orbitalsType <: AbstractArray,
                                orbitalsenergyType <: AbstractArray,
                                occupationType <: AbstractArray}
                            
    niter::Int                                  # Number of iterations
    stopping_criteria::dataType                 # Current stopping criteria
    discretization::discretizationType          # Discretization parameters
    model::modelType                            # Model
    method::methodType                          # Iterative method
    opts::optsType                              # Solver options
    D::densityType                              # Density Matrix at current time
    Dprev::densityType                          # Density Matrix at previous time
    tmpD::densityType                           # Storage for extra Density Matrix (usefull for degeneracy)
    U::orbitalsType                             # Coefficient of orbitals at current time
    ϵ::orbitalsenergyType                       # Orbitals energy at current time
    n::occupationType                           # Occupation number at current time 
    energy::dataType                            # Total energy at current time
    energy_kin::dataType                        # Kinetic energy at current time
    energy_cou::dataType                        # Coulomb energy at current time
    energy_har::dataType                        # Hartree energy at current time
    energy_kin_prev::dataType                   # Kinetic energy at previoustime
    energy_cou_prev::dataType                   # Coulomb energy at previous time
    energy_har_prev::dataType                   # Hartree energy at previous time
    energy_exc::dataType                        # Exchange-correlation energy at current time
    energy_kincor::dataType                     # Kinetic-correlation energy at current time (exists for the LSDA model)
    logbook::logbookType                        # LogBook

    function KhonShamSolver(niter::Int, 
                            stopping_criteria::Real, 
                            discretization::KohnShamDiscretization, 
                            model::AbstractDFTModel, 
                            method::SCFMethod, 
                            opts::SolverOptions, 
                            D::AbstractArray, 
                            Dprev::AbstractArray,
                            tmpD::AbstractArray,
                            U::AbstractArray, 
                            ϵ::AbstractArray, 
                            n::AbstractArray, 
                            energy::Real, 
                            energy_kin::Real, 
                            energy_cou::Real, 
                            energy_har::Real, 
                            energy_kin_prev::Real, 
                            energy_cou_prev::Real, 
                            energy_har_prev::Real,
                            energy_exc::Real, 
                            energy_kincor::Real, 
                            logbook::LogBook)
        new{typeof(discretization),
            typeof(model),
            typeof(method),
            typeof(opts),
            typeof(logbook),
            typeof(energy),
            typeof(D),
            typeof(U),
            typeof(ϵ),
            typeof(n)}(niter, stopping_criteria, discretization, model, method, opts, D, Dprev, tmpD, U, ϵ, n, 
            energy, energy_kin, energy_cou, energy_har, 
            energy_kin_prev, energy_cou_prev, energy_har_prev,
            energy_exc, energy_kincor, logbook)
    end
end