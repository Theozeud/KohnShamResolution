function groundstate(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::AbstractKohnShamResolutionMethod; show_progress = false, kwargs...)
    solver = init(model, discretization, method; kwargs...)
    solve!(solver, method; show_progress = show_progress)
    makesolution(solver)
end

groundstate(problem::DFTProblem; kwargs...) =  groundstate(model(problem), discretization(problem), method(problem); kwargs...)


function init(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::AbstractKohnShamResolutionMethod; lₕ, Nₕ, ε)

    # Check if the method is appropriate
    if !ismethod_for_model(method, model)
        error("This method can't be used for this model.")
    end

    # Init storage array
    D, Dprev    = init_density_matrix(discretization)
    U           = init_coeffs_discretization(discretization)
    ϵ           = init_energy(discretization)
    n           = init_occupation(discretization)
    
    # Init Cache
    cache = init_cache(method, model, discretization; lₕ = lₕ, Nₕ = Nₕ)

    # Registering SolverOptions
    opts = SolverOptions(lₕ,Nₕ,ε)
    


    niter = 0
    current_crit =  2*ε
    val_crit = [current_crit]    

    KhonShamSolver(cache, opts, ϵ, U, n, D, Dprev, niter, val_crit, current_crit)
end

function solve!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod; show_progress = false)
    p = ProgressThresh(solver.opts.ε; enabled = show_progress, desc = "Itération Principale") 
    while solver.current_crit > solver.opts.ε
        println("itérations")
        update!(p, solver.current_crit)
        performstep!(method, solver.cache, solver)
        loopfooter!(solver, method)
    end
end

function loopfooter!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod)
    solver.Rprev .= solver.R
    solver.niter += 1
    solver.current_crit = stopping_criteria(method, solver)
    push!(solver.val_crit, solver.current_crit)
end


function makesolution(solver::KhonShamSolver)
    KohnShamSolution()
end