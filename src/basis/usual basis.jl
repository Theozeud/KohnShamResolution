########################################################################################
#                                   P1 Element
########################################################################################

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

########################################################################################
#                                   P2 Element
########################################################################################

function FunctionP2_node(mesh::OneDMesh, i::Int, T::Type = Float64)
    if i == firstindex(mesh)
        pc = T(mesh[i])
        pr = T(mesh[i+1])
        sr = pr - pc
        right = LaurentPolynomial([T(1)+T(3)*pc/sr + T(2)*(pc/sr)^2, -T(3)/sr - T(4)*pc/sr^2, T(2)/sr^2], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [right], [i], T(0))
    elseif i == lastindex(mesh)
        pl = T(mesh[i-1])
        pc = T(mesh[i])
        sl = pc - pl
        left = LaurentPolynomial([T(1)-T(3)*pc/sl + T(2)*(pc/sl)^2, T(3)/sl - T(4)*pc/sl^2, T(2)/sl^2], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left], [i-1], T(0))
    else
        pl = T(mesh[i-1])
        pc = T(mesh[i])
        pr = T(mesh[i+1])
        sl = pc - pl
        sr = pr - pc
        left = LaurentPolynomial([T(1)-T(3)*pc/sl + T(2)*(pc/sl)^2, T(3)/sl - T(4)*pc/sl^2, T(2)/sl^2], 0, false, T(0))
        right = LaurentPolynomial([T(1)+T(3)*pc/sr + T(2)*(pc/sr)^2, -T(3)/sr - T(4)*pc/sr^2, T(2)/sr^2], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left, right], [i-1,i], T(0))
    end
end

function FunctionP2_mid(mesh::OneDMesh, i::Int, T::Type = Float64)
    @assert i ≠ lastindex(mesh) "The mid-point function can't be defined on the mid-point right after the last point of the mesh."
    pl = T(mesh[i])
    pr = T(mesh[i+1])
    pc = (pr + pl)/2
    s2 = (pr - pl)^2
    p = LaurentPolynomial([1-4*pc^2/s2, 8*pc/s2, -4/s2], 0, false, T(0))
    PiecewiseLaurentPolynomial(mesh, [p], [i], T(0))
end

function P2Basis(mesh::OneDMesh{TM}, T::Type = Float64; left::Bool = true, right::Bool = left) where TM
    basis = PiecewiseLaurentPolynomial{T,TM}[]
    if left
        push!(basis, FunctionP2_node(mesh, firstindex(mesh), T))
    end
    for i ∈ eachindex(mesh)[begin+1:end-1]
        push!(basis, FunctionP2_mid(mesh, i-1, T))
        push!(basis, FunctionP2_node(mesh, i, T))
    end
    push!(basis, FunctionP2_mid(mesh, lastindex(mesh) - 1, T))
    if right
        push!(basis, FunctionP2_node(mesh, lastindex(mesh), T))
    end
    LaurentPolynomialBasis(basis)
end