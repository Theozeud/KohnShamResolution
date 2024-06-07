function HatFunctionP1(mesh::OneDMesh, i::Int, T::Type = Float64)
    if i == firstindex(mesh)
        pc = mesh[i]
        pr = mesh[i+1]
        right = LaurentPolynomial([pr/(pr-pc),T(1)/(pc-pr)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [right], [i], T(0))
    elseif i == lastindex(mesh)
        pl = mesh[i-1]
        pc = mesh[i]
        left = LaurentPolynomial([pl/(pl-pc),T(1)/(pc-pl)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left], [i-1], T(0))
    else
        pl = mesh[i-1]
        pc = mesh[i]
        pr = mesh[i+1]
        left = LaurentPolynomial([pl/(pl-pc),T(1)/(pc-pl)], 0, false, T(0))
        right = LaurentPolynomial([pr/(pr-pc),T(1)/(pc-pr)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left, right], [i-1,i], T(0))
    end
end


function HatBasis(mesh::OneDMesh, T::Type = Float64)
    
end