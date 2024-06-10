abstract type AbstractKohnShamSolution end

struct KohnShamSolution
    success::String 
    ϵ
    n
    niter::Int
    crit::Vector

    function KohnShamSolution(solver::KhonShamSolver)
        if solver.niter == solver.opts.maxiter
            success = "MAXITERS"
        else
            success = "SUCCESS"
        end
        new(success, solver.ϵ[] , solver.n, solver.niter, solver.values_stop_crit)
        
    end
end

function Base.show(io::IO, sol::KohnShamSolution)
    printstyled(io, "Sucess = "; bold = true)
    if sol.success == "SUCCESS"
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :green)
    else
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :red)
    end
    printstyled(io, "ϵ = "; bold = true)
    println(io, string(sol.ϵ))
    printstyled(io, "n = "; bold = true)
    println(io, string(sol.n))
    printstyled(io, "niter = "; bold = true)
    println(io, string(sol.niter))
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, string(last(sol.crit)))
end

# afficher en plus le fondamentale et l'orbital correspondant, 
# affciher les énergies ϵ seulement pour les n non nuls
# retenir les états occupés