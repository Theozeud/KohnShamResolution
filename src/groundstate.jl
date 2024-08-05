function groundstate(model::AbstractDFTModel, discretization::KohnShamDiscretization, method::AbstractKohnShamResolutionMethod; show_progress = false, kwargs...)
    solver = init(model, discretization, method; kwargs...)
    solve!(solver; show_progress = show_progress)
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
    degen_tol::Real = eps(bottom_type(discretization.basis)),
    potential = :pde)

    # Valid option for resol_hartree
    @assert potential ∈ [:pde, :integral] "potential must be :pde or :integral"

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
    
    #  SolverOptions
    opts = SolverOptions(tol, maxiter, quad_method, T(quad_reltol), T(quad_abstol), T(hartree), degen_tol, potential)
    niter = 0
    current_stop_crit =  2*T(tol)
    values_stop_crit = T[]    
    Ehisto = T[]

    KhonShamSolver(discretization, model, method, D, Dprev, U, ϵ, n, Ehisto, niter, values_stop_crit, current_stop_crit, opts)
end

function solve!(solver::KhonShamSolver; show_progress = false)
    p = ProgressThresh(solver.opts.ε; enabled = show_progress, desc = "Itération Principale") 
    while solver.current_stop_crit > solver.opts.ε && solver.niter < solver.opts.maxiter
        update!(p, solver.current_stop_crit)
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

    # Reset Solver
    reset_cache!(solver)
end


function loopfooter!(solver::KhonShamSolver)
    push!(solver.Ehisto, min(solver.ϵ...))
    solver.current_stop_crit = stopping_criteria(solver)
    push!(solver.values_stop_crit, solver.current_stop_crit)
    solver.Dprev  .= solver.D
    solver.niter += 1
end

#function stopping_criteria(solver::KhonShamSolver)
#    m = solver.discretization.mesh
#    sqrt((m[end] - m[begin])*sum([abs(solver.D(x)- solver.Dprev(x))^2 for x ∈ m[begin:end-1]])/(length(m)-1))
#end
function stopping_criteria(solver::KhonShamSolver)
    norm(solver.D - solver.Dprev)
end
    
function makesolution(solver::KhonShamSolver)
    KohnShamSolution(solver)
end

function reset_cache!(solver::KhonShamSolver)
    @unpack tmp_D, tmp_Dstar, tmp_U, Hartree, tmp_exc, tmp_ϵ, tmp_n = solver.discretization.cache
    solver.D .= tmp_D
    solver.U .= tmp_U
    solver.ϵ .= tmp_ϵ
    solver.n .= tmp_n 
    tmp_D       .= zero(tmp_D)
    tmp_Dstar   .= zero(tmp_Dstar)
    tmp_U       .= zero(tmp_U)
    tmp_exc     .= zero(tmp_exc)
    Hartree     .= zero(Hartree)
    tmp_ϵ       .= zero(tmp_ϵ)
    tmp_n       .= zero(tmp_n)
end