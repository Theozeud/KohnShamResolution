# Structure Solution 

const L_QUANTUM_LABELS = ("s", "p", "d", "f", "g", "h", "i")
const SPIN_LABELS = ("↑", "↓")

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
    Energy_kincor           # Final Kinetic-correlation energy

    occupation_number       # Final occupation number
    orbitals_energy         # Final orbitals energy
    orbitals                # Final coefficients of orbitals
    density_coeffs          # Final coefficients of density
    
    log::LogBook            # LogBook of all recorded quantities

    name::String            # Name of the solution

    function KohnShamSolution(solver::KhonShamSolver, name::String)

        problem = DFTProblem(solver.model, solver.discretization, solver.method)

        success = solver.niter == solver.opts.maxiter ? "MAXITERS" : "SUCCESS"

        occupation_number = make_occupation_number(solver.discretization, solver)

        new(problem, success, solver.opts, solver.niter, solver.stopping_criteria, 
            solver.energy,solver.energy_kin, solver.energy_cou, solver.energy_har, solver.energy_exc, solver.energy_kincor,
            occupation_number, solver.ϵ, solver.U, solver.D,solver.logbook, name)
    end
end


# Make occupation number

function make_occupation_number(::LDADiscretization, solver::KhonShamSolver)
    @unpack ϵ, n = solver
    index = findall(x->x ≠ 0, n)
    index_sort = sortperm(solver.ϵ[index])
    new_index = index[index_sort]
    return [(string(i[2]+ i[1] -1, L_QUANTUM_LABELS[i[1]]), solver.ϵ[i], solver.n[i]) for i ∈ new_index]
end

function make_occupation_number(::LSDADiscretization, solver::KhonShamSolver)
    @unpack ϵ, n = solver
    index = findall(x->x ≠ 0, n)
    index_sort = sortperm(solver.ϵ[index])
    new_index = index[index_sort]
    return [(string(i[2]+ i[1] -1, L_QUANTUM_LABELS[i[1]],SPIN_LABELS[i[3]]), solver.ϵ[i], solver.n[i]) for i ∈ new_index]
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
        display_occupation_number(io, sol.problem.discretization, sol.occupation_number[i])
    end
    printstyled(io, "niter = "; bold = true)
    println(io, string(sol.niter))
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, string(sol.stopping_criteria ))
end

function display_occupation_number(io::IO, ::LDADiscretization, occupation_number)
    printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
end

function display_occupation_number(io::IO, ::LSDADiscretization, occupation_number)
    printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
end

# Compute eigenvector
function eigenvector(sol::KohnShamSolution, n::Int, l::Int, σ::Int, x)
    eigenvector(sol.problem.discretization, sol, n, l, σ, x)
end

function eigenvector(discretization::LDADiscretization, sol::KohnShamSolution, n::Int, l::Int, σ::Int, x)
    @assert 0 ≤ l ≤ n-1
    tmp = discretization.basis(sol.orbitals[l+1,:, n-l], x)
    if iszero(x) && tmp ≈ zero(tmp)
        return zero(tmp)
    else
        return  1/sqrt(4π * x) * tmp 
    end
end

function eigenvector(discretization::LSDADiscretization, sol::KohnShamSolution, n::Int, l::Int, σ::Int, x)
    @assert 0 ≤ l ≤ n-1
    tmp = discretization.basis(sol.orbitals[l+1,:, n-l,σ], x)
    if iszero(x) && tmp ≈ zero(tmp)
        return zero(tmp)
    else
        return  1/sqrt(4π * x) * tmp 
    end
end


# Compute density
function density(sol::KohnShamSolution, x::Real)
    compute_density(sol.problem.discretization, sol.density_coeffs, x)
end

function density(sol::KohnShamSolution, X::AbstractVector)
    [compute_density(sol.problem.discretization, sol.density_coeffs, x) for x ∈ X]
end

# Check integral of density
function total_charge(sol::KohnShamSolution)
    @unpack Rmax = sol.problem.discretization
    f(x,p) = density(sol, x) * x^2
    prob = IntegralProblem(f, (zero(Rmax),Rmax))
    return 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
end