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
    eigen(inv(B) * A)
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
function approximate_integral(f, domain; method = QuadGKJL(), reltol = 1e-3, abstol = 1e-3)
    prob = IntegralProblem((x,p) -> f(x), domain)
    solve(prob, method; reltol = reltol, abstol = abstol).u
end