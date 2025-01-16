# Structure Solution 

const L_QUANTUM_LABELS = ("s", "p", "d", "f", "g", "h", "i")

struct KohnShamSolution
    problem                 # Problem solved

    success::String         # Print the final state of the solver
                            #        Can be : SUCCESS or MAXITERS

    solveropts              # Option of the solver used
    
    niter::Int              # Number of iterations
    stopping_criteria       # Final stopping criteria
              
    Energy                  # Final energy
    Energy_kin              # Final Kinertic energy
    Energy_cou              # Final Coulomb energy
    Energy_har              # Final Hartree energy
    Energy_exc              # Final Exchange-correlation energy

    occupation_number       # Final occupation number
    orbitals_energy         # Final orbitals energy
    orbitals                # Final coefficients of orbitals
    density_coeffs          # Final coefficients of density
    
    log::LogBook            # LogBook of all recorded quantities

    name::String            # Name of the solution

    function KohnShamSolution(solver::KhonShamSolver, name::String)

        problem = DFTProblem(solver.model, solver.discretization, solver.method)

        success = solver.niter == solver.opts.maxiter ? "MAXITERS" : "SUCCESS"

        @unpack ϵ, n = solver

        index = findall(x->x ≠ 0, n)
        ϵ_full = solver.ϵ[index]
        index_sort = sortperm(ϵ_full)
        new_index = index[index_sort]
        occupation_number = [(string(i[2]+ i[1] -1, L_QUANTUM_LABELS[i[1]]), solver.ϵ[i], 
                                solver.n[i]) for i ∈ new_index]

        new(problem, success, solver.opts, solver.niter, solver.stopping_criteria, 
            solver.energy,solver.energy_kin, solver.energy_cou, solver.energy_har, solver.energy_exc, 
            occupation_number, ϵ, solver.U, solver.D, solver.logbook, name)
    end
end

# Show function to print a summary of the solution in the stream

function Base.show(io::IO, sol::KohnShamSolution)
    printstyled(io, "Name : " * (sol.name) * "\n"; bold = true)
    printstyled(io, "Sucess = "; bold = true)
    if sol.success == "SUCCESS"
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :green)
    else
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :red)
    end
    printstyled(io, "Energy = $(sol.Energy) \n"; bold = true, color = :green)
    printstyled(io, "Occupation number = \n"; bold = true, color = :blue)
    for i ∈ eachindex(sol.occupation_number)
        printstyled(io, "            $(sol.occupation_number[i][1]) : ($(sol.occupation_number[i][2]),$(sol.occupation_number[i][3])) \n"; bold = true, color = :blue)
    end
    printstyled(io, "niter = "; bold = true)
    println(io, string(sol.niter))
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, string(sol.stopping_criteria ))
end


# Compute eigenvector
function eigenvector(sol::KohnShamSolution, n::Int, l::Int, x)
    @assert 0 ≤ l ≤ n-1
    tmp = sol.problem.discretization.basis(sol.orbitals[l+1,:, n-l], x)
    if iszero(x) && tmp ≈ zero(tmp)
        return zero(tmp)
    else
        return  1/sqrt(4π * x) * tmp 
    end
end

# Compute density
function density(sol::KohnShamSolution, x)
    compute_density(sol.problem.discretization, sol.density_coeffs, x)
end



