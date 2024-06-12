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
    hartree::Bool = true)

    # Init storage array
    D, Dprev    = init_density_matrix(discretization)
    U           = init_coeffs_discretization(discretization)
    ϵ           = init_energy(discretization)
    n           = init_occupation(discretization)
    
    # Init Cache
    cache = init_cache(method, model, discretization)

    # Registering SolverOptions
    opts = SolverOptions(tol, maxiter, quad_method, quad_reltol, quad_abstol, hartree)
    niter = 0
    current_stop_crit =  2*tol
    values_stop_crit = [current_stop_crit]    

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