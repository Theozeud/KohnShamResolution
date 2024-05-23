# Quadrature


# Solving linear problem
```
    Solve Linear Problem of the form AX = b
    where X is the unknown
```
function solve_linear_problem(A::AbstractMatrix, b::AbstractVector)
    A\b
end


# Solve Generalized Eigenvalue problem


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
