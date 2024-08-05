# Solving linear problem
function solve_linear_problem(A::AbstractMatrix, b::AbstractVector)
    A\b
end

# Solve Generalized Eigenvalue problem
function solve_generalized_eigenvalue_problem(A,B)
    eigen(inv(B) * A)
end

# Approximate integral
function approximate_integral(f, domain; method = QuadGKJL(), reltol = 1e-3, abstol = 1e-3)
    prob = IntegralProblem((x,p) -> f(x), domain)
    solve(prob, method; reltol = reltol, abstol = abstol).u
end
