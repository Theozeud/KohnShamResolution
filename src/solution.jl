abstract type AbstractKohnShamSolution end

struct KohnShamSolution
    success::String 
    niter::Int
    crit::Vector
    Energy::Real
    occup::Vector
    ϵ
    ϵhisto
    Energyhisto
    ρ

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
        ρ = solver.opts.light ? nothing : build_density2!(solver.discretization, solver.D)

        new(success, solver.niter, solver.values_stop_crit, last(solver.Energyhisto), occup, ϵ, solver.ϵhisto, solver.Energyhisto, ρ)
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
    printstyled(io, "Energy = "; bold = true, color = :green)
    printstyled(io, sol.Energy; bold = true, color = :green)
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