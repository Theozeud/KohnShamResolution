
function groundstate(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::SCFMethod; name = "", kwargs...)
    solver = init(model, discretization, method; kwargs...)
    solve!(solver)
    makesolution(solver, name)
end

groundstate(problem::DFTProblem; kwargs...) =  groundstate(model(problem), discretization(problem), method(problem); kwargs...)

function init(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::SCFMethod; 
    scftol::Real, 
    maxiter::Int = 100,
    quad_method = QuadGKJL(),
    quad_reltol::Real  = 1e-3,
    quad_abstol::Real   = 1e-3,
    hartree::Real = 1, 
    degen_tol::Real = eps(bottom_type(discretization.basis)),
    logconfig = LogConfig())

    # Set the type of number as the one of the discretization basis
    T = discretization.elT

    # Init Cache of the Discretisation
    init_cache!(discretization, model, hartree)

    # Init storage array
    D               = init_density_matrix(discretization)
    Dprev           = init_density_matrix(discretization)
    U               = init_coeffs_discretization(discretization)
    ϵ               = init_energy(discretization)
    n               = init_occupation(discretization)
    energy          = zero(T)
    energy_kin      = zero(T)
    energy_cou      = zero(T)
    energy_har      = zero(T)
    energy_kin_prev = zero(T)
    energy_cou_prev = zero(T)
    energy_har_prev = zero(T)
    energy_exc      = zero(T)
    energy_kincor   = zero(T)
    
    #  SolverOptions
    opts = SolverOptions(T(scftol), maxiter, quad_method, T(quad_reltol), T(quad_abstol), T(hartree), T(degen_tol))

    # Init log parameters
    niter = 0
    stopping_criteria = zero(T)
    logbook = LogBook(logconfig, T)
    
    KhonShamSolver(niter, stopping_criteria, discretization, model, method, opts, D, Dprev, U, ϵ, n, 
                   energy, energy_kin, energy_cou, energy_har, 
                   energy_kin_prev, energy_cou_prev, energy_har_prev,
                   energy_exc, energy_kincor, logbook)
end


function solve!(solver::KhonShamSolver)
    while (solver.stopping_criteria > solver.opts.scftol || iszero(solver.niter)) && solver.niter < solver.opts.maxiter
        println("Iteration : $(solver.niter)")
        loopheader!(solver)
        performstep!(solver)
        loopfooter!(solver)
    end
    compute_energy!(solver.discretization, solver)
end

function performstep!(solver::KhonShamSolver)
    # STEP 1 : Resolution of the generalized eigenvalue problem to find atomic orbitals and corresonding energies
    find_orbital!(solver.discretization, solver)

    # STEP 2 : Build the n matrix using the Aufbau principle
    aufbau!(solver)

    # STEP 3 : Normaization of eigenvectors
    # This is done after aufbau to normalize only eigenvectors we need
    normalization!(solver.discretization, solver)

    # STEP 4 : Compute density star
    density_matrix!(solver.discretization, solver)

    # STEP 5 : Compute energy
    compute_energy!(solver.discretization, solver)

    # STEP 6 : Compute new density
    update_density!(solver.method, solver)

    # STEP 7 : Update Solver
    update_solver!(solver)
end

function loopheader!(solver::KhonShamSolver)
    solver.Dprev            .= solver.D
    solver.D                .= zero(solver.D)
    solver.n                .= zero(solver.n)
    solver.energy_kin_prev   = solver.energy_kin
    solver.energy_cou_prev   = solver.energy_cou
    solver.energy_har_prev   = solver.energy_har                                          
end 

function loopfooter!(solver::KhonShamSolver)
    solver.stopping_criteria = stopping_criteria(solver)            # COMPUTE THE NEW STOPPING CRITERIA                                 
    solver.niter += 1                                               # INCREASE THE NUMBER OF ITERATIONS DONE
    update_log!(solver)                                             # UPDATE THE LOG
end 


function stopping_criteria(solver::KhonShamSolver)
    norm(solver.D - solver.Dprev)
end
    
function makesolution(solver::KhonShamSolver, name::String)
    KohnShamSolution(solver, name)
end

function update_solver!(solver::KhonShamSolver)
    @unpack Cprev, C, Cᵨ, Cᵨprev = solver.discretization.cache
    solver.discretization.cache.Cprev = solver.discretization.cache.C
    solver.discretization.cache.Cᵨprev = Cᵨ
    nothing
end

function update_log!(solver::KhonShamSolver)
    @unpack logbook = solver
    @unpack occupation_number, orbitals_energy, stopping_criteria, energy, density = logbook.config

    # STORE THE STOPPING CRITERIA 
    stopping_criteria ? push!(solver.logbook.stopping_criteria_log, solver.stopping_criteria)  : nothing     
    
    # STORE THE ORBITALS ENERGY
    orbitals_energy ?   push!(solver.logbook.orbitals_energy_log, copy(solver.ϵ))              : nothing     
    
    # STORE THE TOTAL ENERGY
    if energy
        if isthereExchangeCorrelation(solver.model)
            compute_exchangecorrelation_energy!(solver.discretization, solver)
        end
        compute_total_energy!(solver.discretization, solver)
        push!(solver.logbook.energy_log, solver.energy)                           
    end
end