struct SolverOptions{T}
    scftol::T                       # SCF tolerance
    maxiter::Int                    # Maximum of iteration done

    quad_method                     # Method to use for quadrature of integrals
    quad_reltol::T                  # Relativ tolerance for the quadrature of integrals
    quad_abstol::T                  # Absolute tolerance for the quadrature of integrals

    hartree::T                      # Coefficient multiply to the Hartree Matrix : 
                                    # 0 -> no hartree term, 1-> full hartree term

    degen_tol::T                    # Tolerance to consider degenescence of orbitals energy

    verbose::UInt8                  # Say how many details are printed at the end 
                                    # of each iterations :
                                    # 0 : Zero details
                                    # 1 : Iterations
                                    # 2 : Computations                
end


mutable struct KohnShamSolver{  discretizationType <: KohnShamDiscretization,
                                modelType <: AbstractDFTModel,
                                methodType <: SCFMethod,
                                cacheType <: SCFCache,
                                optsType <: SolverOptions,
                                logbookType <: LogBook,
                                dataType <: Real}
                            
    niter::Int                                  # Number of iterations
    stopping_criteria::dataType                 # Current stopping criteria
    discretization::discretizationType          # Discretization parameters
    model::modelType                            # Model
    method::methodType                          # Iterative method
    cache::cacheType                            # Cache associated to the method
    opts::optsType                              # Solver options
    energies::Dict{Symbol,dataType}             # Storages of the energies
    logbook::logbookType                        # LogBook

    #=
    energy::dataType                            # Total energy at current time
    energy_kin::dataType                        # Kinetic energy at current time
    energy_cou::dataType                        # Coulomb energy at current time
    energy_har::dataType                        # Hartree energy at current time
    energy_exc::dataType                        # Exchange-correlation energy at current time
    energy_kincor::dataType                     # Kinetic-correlation energy at current time (exists for the LSDA model)
    energy_prev::dataType                       # Total energy at previous time
    energy_kin_prev::dataType                   # Kinetic energy at previoustime
    energy_cou_prev::dataType                   # Coulomb energy at previous time
    energy_har_prev::dataType                   # Hartree energy at previous time
    =#
    
    
    function KohnShamSolver(niter::Int, 
                            stopping_criteria::Real, 
                            discretization::KohnShamDiscretization, 
                            model::AbstractDFTModel, 
                            method::SCFMethod,
                            cache::SCFCache,
                            opts::SolverOptions,  
                            energies::Dict{Symbol,<:Real}, 
                            logbook::LogBook)
        new{typeof(discretization),
            typeof(model),
            typeof(method),
            typeof(cache),
            typeof(opts),
            typeof(logbook),
            typeof(stopping_criteria)}( niter, 
                                        stopping_criteria, 
                                        discretization, 
                                        model, 
                                        method, 
                                        cache, 
                                        opts,
                                        energies,
                                        logbook)
    end
end