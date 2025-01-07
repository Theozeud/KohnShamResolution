
function groundstate(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::SCFMethod; kwargs...)
    solver = init(model, discretization, method; kwargs...)
    solve!(solver)
    makesolution(solver)
end

groundstate(problem::DFTProblem; kwargs...) =  groundstate(model(problem), discretization(problem), method(problem); kwargs...)

function init(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::SCFMethod; 
    scftol::Real, 
    maxiter::Int = 100,
    quad_method = QuadGKJL(),
    quad_reltol::Real  = 1e-3,
    quad_abstol::Real   = 1e-3,
    hartree::Real = 0, 
    degen_tol::Real = eps(bottom_type(discretization.basis)),
    logconfig = LogConfig())

    # Set the type of number as the one of the discretization basis
    T = discretization.elT

    # Init Cache of the Discretisation
    init_cache!(discretization, model, hartree)

    # Init storage array
    D           = init_density_matrix(discretization)
    Dprev       = init_density_matrix(discretization)
    U           = init_coeffs_discretization(discretization)
    ϵ           = init_energy(discretization)
    n           = init_occupation(discretization)
    energy      = zero(T)
    
    #  SolverOptions
    opts = SolverOptions(T(scftol), maxiter, quad_method, T(quad_reltol), T(quad_abstol), T(hartree), T(degen_tol))

    # Init log parameters
    niter = 0
    stopping_criteria = zero(T)
    logbook = LogBook(logconfig, T)
    
    KhonShamSolver(niter, stopping_criteria, discretization, model, method, opts, D, Dprev, U, ϵ, n, energy, logbook)
end



function solve!(solver::KhonShamSolver)
    while (solver.stopping_criteria > solver.opts.scftol || iszero(solver.niter)) && solver.niter < solver.opts.maxiter
        performstep!(solver)
        loopfooter!(solver)
    end
end

function performstep!(solver::KhonShamSolver)

    # STEP 1 : Resolution of the generalized eigenvalue problem to find atomic orbitals and corresonding energies
    find_orbital!(solver.discretization, solver)

    # STEP 2 : Build the n matrix using the Aufbau principle
    aufbau!(solver)

    # STEP 3 : Compute density star
    density_matrix!(solver.discretization)

    # STEP 4 : Compute new density
    update_density!(solver.method, solver)

    # STEP 5 : Update Solver
    update_solver!(solver)
end


function loopfooter!(solver::KhonShamSolver)
    solver.stopping_criteria = stopping_criteria(solver)            # COMPUTE THE NEW STOPPING CRITERIA 
    solver.Dprev  .= solver.D                                       # UPDATE THE CURRENT DENSITY
    solver.niter += 1                                               # INCREASE THE NUMBER OF ITERATIONS DONE
    update_log!(solver)                                             # UPDATE THE LOG
end 


function stopping_criteria(solver::KhonShamSolver)
    norm(solver.D - solver.Dprev)
end
    
function makesolution(solver::KhonShamSolver)
    KohnShamSolution(solver)
end

function update_solver!(solver::KhonShamSolver)
    @unpack tmp_D, tmp_Dstar, tmp_U, tmp_ϵ, tmp_n = solver.discretization.tmp_cache
    solver.D .= tmp_D
    solver.U .= tmp_U
    solver.ϵ .= tmp_ϵ
    solver.n .= tmp_n 
    tmp_D       .= zero(tmp_D)
    tmp_Dstar   .= zero(tmp_Dstar)
    nothing
end

function update_log!(solver::KhonShamSolver)
    @unpack logbook = solver
    @unpack occupation_number, orbitals_energy, stopping_criteria, energy, density = logbook.config
    stopping_criteria ? push!(solver.logbook.stopping_criteria_log, solver.stopping_criteria)  : nothing     # STORE THE STOPPING CRITERIA 
    orbitals_energy ?   push!(solver.logbook.orbitals_energy_log, copy(solver.ϵ[1,:]))         : nothing     # STORE THE ORBITALS ENERGY
    energy ?            push!(solver.logbook.energy_log, solver.discretization.cache.Energy)   : nothing     # STORE THE TOTAL ENERGY
    nothing
end