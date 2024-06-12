function HatFunctionP1(mesh::OneDMesh, i::Int, T::Type = Float64)
    if i == firstindex(mesh)
        pc = T(mesh[i])
        pr = T(mesh[i+1])
        right = LaurentPolynomial([pr/(pr-pc),T(1)/(pc-pr)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [right], [i], T(0))
    elseif i == lastindex(mesh)
        pl = T(mesh[i-1])
        pc = T(mesh[i])
        left = LaurentPolynomial([pl/(pl-pc),T(1)/(pc-pl)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left], [i-1], T(0))
    else
        pl = T(mesh[i-1])
        pc = T(mesh[i])
        pr = T(mesh[i+1])
        left = LaurentPolynomial([pl/(pl-pc),T(1)/(pc-pl)], 0, false, T(0))
        right = LaurentPolynomial([pr/(pr-pc),T(1)/(pc-pr)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left, right], [i-1,i], T(0))
    end
end


function HatBasis(mesh::OneDMesh, T::Type = Float64; left::Bool = true, right::Bool = left)
    if left && right
        index = eachindex(mesh)
    elseif !left && !right
        index = eachindex(mesh)[begin+1:end-1]
    elseif !left
        index = eachindex(mesh)[begin+1:end]
    else
        index = eachindex(mesh)[begin:end-1]
    end
    LaurentPolynomialBasis([HatFunctionP1(mesh, i, T) for i ∈ index])
end