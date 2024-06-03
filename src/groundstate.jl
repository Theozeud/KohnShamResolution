function groundstate(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; show_progress = false, kwargs...)
    solver = init(model, method; kwargs...)
    solve!(solver, method; show_progress = show_progress)
    makesolution(solver)
end

groundstate(problem::DFTProblem; kwargs...) =  groundstate(model(problem), method(problem); kwargs...)


function init(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; lₕ, Nₕ, ε)

    # Check if the method is appropriate
    if !ismethod_for_model(method, model)
        error("This method can't be used for this model.")
    end

    # Init Cache
    cache = init_cache(model, method; lₕ = lₕ, Nₕ = Nₕ)

    # Registering SolverOptions
    opts = SolverOptions(lₕ,Nₕ,ε)
    
    # Init storage array
    ϵ       = zeros(lₕ+1,(2lₕ+1)Nₕ)
    U       = zeros(lₕ+1, (2lₕ+1)Nₕ, Nₕ)
    n       = zeros(lₕ+1, (2lₕ+1)Nₕ)
    R       = zeros(lₕ+1, Nₕ, Nₕ)
    Rprev   = zeros(lₕ+1, Nₕ, Nₕ)

    niter = 0
    current_crit =  stopping_criteria(method, R, Rprev)
    val_crit = [current_crit]

    solver(cache, opts, ϵ, U, n, R, Rprev, niter, val_crit, current_crit)
end

function solve!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod; show_progress = false)
    p = ProgressThresh(solver.opts.ϵ; enabled = show_progress, desc = "Itération Principale")
    while solver.current_crit > solver.opts.ϵ
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