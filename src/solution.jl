abstract type AbstractKohnShamSolution end

struct KohnShamSolution
    success::String 
    occup::Vector
    E::Real
    ϵ::Vector
    eigvects
    niter::Int
    crit::Vector

    function KohnShamSolution(solver::KhonShamSolver)

        if solver.niter == solver.opts.maxiter
            success = "MAXITERS"
        else
            success = "SUCCESS"
        end
        @unpack ϵ, U, n = solver

        index = findall(x->x ≠ 0, n)
        ϵ_full = ϵ[index]
        index_sort = sortperm(ϵ_full)
        new_index = index[index_sort]
        
        occup = [(string(i[2], convert_into_l(i[1]-1)), ϵ[i], n[i]) for i ∈ new_index]
        eigvects = build_eigenvector(solver.discretization, solver.U; Index =  new_index)

        new(success, occup, ϵ[first(new_index)], sort(vec(ϵ)), eigvects, solver.niter, solver.values_stop_crit)
    end
end

function Base.show(io::IO, sol::KohnShamSolution)
    printstyled(io, "Sucess = "; bold = true)
    if sol.success == "SUCCESS"
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :green)
    else
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :red)
    end
    printstyled(io, "ϵn = "; bold = true)
    println(io, string(sol.occup))
    printstyled(io, "E = "; bold = true, color = :green)
    printstyled(io, sol.E; bold = true, color = :green)
    println("")
    printstyled(io, "niter = "; bold = true)
    println(io, string(sol.niter))
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, string(last(sol.crit)))
end

function convert_into_l(l::Int)
    if l == 0
        return "s"
    elseif l == 1
        return "p"
    elseif l == 2
        return "d"
    elseif l == 3
        return "f" 
    elseif l == 4
        return "g" 
    elseif l == 5
        return "h" 
    elseif l == 6
        return "i"
    else
        return "l_"*string(l)
    end
end 

Base.length(sol::KohnShamSolution) = length(sol.occup)
Base.eachindex(sol::KohnShamSolution) = eachindex(sol.occup)

eigenvalue(sol::KohnShamSolution, i::Int) = sol.occup[i][2]
eigenvector(sol::KohnShamSolution, i::Int) = sol.eigvects[i]
occup(sol::KohnShamSolution, i::Int) = sol.occup[i][3]

eigenvalue(sol::KohnShamSolution) = [sol.occup[i][2] for i in eachindex(sol)]
eigenvector(sol::KohnShamSolution) = [sol.eigvects[i] for i in eachindex(sol)]
occup(sol::KohnShamSolution) = [sol.occup[i][3] for i in eachindex(sol)]

# afficher en plus le fondamentale et l'orbital correspondant, 
# affciher les énergies ϵ seulement pour les n non nuls
# retenir les états occupés