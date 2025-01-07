
struct KohnShamSolution
    success::String         # Print the final state of the solver
                            #        Can be : SUCCESS or MAXITERS
    
    niter::Int              # Number of iterations
    stopping_criteria       # Final stopping criteria
              
    Energy                  # Final Energy
    occupation_number       # Final occupation number
    orbitals_energy         # Final orbitals energy
    
    log::LogBook            # LogBook of all recorded quantities

    function KohnShamSolution(solver::KhonShamSolver)

        success = solver.niter == solver.opts.maxiter ? "MAXITERS" : "SUCCESS"

        @unpack ϵ, n = solver

        index = findall(x->x ≠ 0, n)
        ϵ_full = solver.ϵ[index]
        index_sort = sortperm(ϵ_full)
        new_index = index[index_sort]
        occupation_number = [(string(i[2], convert_into_l(i[1]-1)), solver.ϵ[i], solver.n[i]) for i ∈ new_index]

        new(success, solver.niter, solver.stopping_criteria, solver.energy, occupation_number, ϵ, solver.logbook)
    end
end

function Base.show(io::IO, sol::KohnShamSolution)
    printstyled(io, "Sucess = "; bold = true)
    if sol.success == "SUCCESS"
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :green)
    else
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :red)
    end
    printstyled(io, "Energy = "; bold = true, color = :green)
    printstyled(io, sol.Energy; bold = true, color = :green)
    println("")
    printstyled(io, "niter = "; bold = true)
    println(io, string(sol.niter))
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, string(sol.stopping_criteria ))
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