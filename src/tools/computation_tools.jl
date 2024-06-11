# Solving linear problem
```
    Solve Linear Problem of the form AX = b
    where X is the unknown
```
function solve_linear_problem(A::AbstractMatrix, b::AbstractVector)
    A\b
end


# Solve Generalized Eigenvalue problem
function solve_generalized_eigenvalue_problem(A,B)
    eigen(A,B)
end

# Optimisation


# Tensor product
function tensorproduct(X::AbstractVector, Y::AbstractVector)
    XY = zeros(size(X)[1],size(Y)[1])
    for i ∈ eachindex(X)
        for j ∈ eachindex(Y)
            XY[i,j] = X[i]Y[j]
        end 
    end
    XY
end


# Approximate integral
function approximate_integral(f, domain; method = QuadGKJL(), retol = 1e-3, abstol = 1e-3)
    prob = IntegralProblem((x,p) -> f, domain, p)
    solve(prob, method, retol = retol, abstol = abstol).u
end