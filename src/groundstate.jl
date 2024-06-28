function groundstate(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::AbstractKohnShamResolutionMethod; show_progress = false, kwargs...)
    solver = init(model, discretization, method; kwargs...)
    solve!(solver, method; show_progress = show_progress)
    makesolution(solver)
end

groundstate(problem::DFTProblem; kwargs...) =  groundstate(model(problem), discretization(problem), method(problem); kwargs...)

function init(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::AbstractKohnShamResolutionMethod; 
    tol::Real, 
    maxiter::Int = 100,
    quad_method = QuadGKJL(),
    quad_reltol::Real  = 1e-3,
    quad_abstol::Real   = 1e-3,
    hartree::Real = 0, 
    degen_tol::Real = eps(bottom_type(discretization.basis)))

    # Set the type of number as the one of the discretization basis
    T = bottom_type(discretization.basis)

    # Init storage array
    D, Dprev    = init_density_matrix(discretization, T)
    U           = init_coeffs_discretization(discretization, T)
    ϵ           = init_energy(discretization, T)
    n           = init_occupation(discretization, T)
    
    # Init Cache
    cache = init_cache(method, model, discretization)

    # Registering SolverOptions
    opts = SolverOptions(tol, maxiter, quad_method, T(quad_reltol), T(quad_abstol), T(hartree), degen_tol)
    niter = 0
    current_stop_crit =  2*T(tol)
    values_stop_crit = T[]    

    KhonShamSolver(discretization, model, D, Dprev, U, ϵ, n, niter, values_stop_crit, current_stop_crit, cache, opts)
end

function solve!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod; show_progress = false)
    p = ProgressThresh(solver.opts.ε; enabled = show_progress, desc = "Itération Principale") 
    while solver.current_stop_crit > solver.opts.ε && solver.niter < solver.opts.maxiter
        update!(p, solver.current_stop_crit)
        performstep!(method, solver)
        loopfooter!(solver, method)
    end
end

function loopfooter!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod)
    solver.current_stop_crit = stopping_criteria(method, solver)
    push!(solver.values_stop_crit, solver.current_stop_crit)
    solver.Dprev  = solver.D
    solver.niter += 1
end

function makesolution(solver::KhonShamSolver)
    KohnShamSolution(solver)
end