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

P1Basis(mesh::OneDMesh, T::Type = Float64; left::Bool = true, right::Bool = left) = HatBasis(mesh, T; left = left, right = right)

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

########################################################################################
#                             Peano hierarchical shape function
########################################################################################

@memoize function Bubble_laurent(mesh::OneDMesh, i::Int, j::Int, T::Type = Float64)
    @assert i ≤ lastindex(mesh) - 1
    pl = T(mesh[i])
    pr = T(mesh[i+1])
    if j == 1
        return LaurentPolynomial([pr/(pr-pl),T(1)/(pl-pr)], 0, false, T(0))
    elseif j == 2
        return LaurentPolynomial([pl/(pl-pr),T(1)/(pr-pl)], 0, false, T(0))
    elseif j == 3
        p₁ = LaurentPolynomial([pr/(pr-pl),T(1)/(pl-pr)], 0, false, T(0))
        p₂ = LaurentPolynomial([pl/(pl-pr),T(1)/(pr-pl)], 0, false, T(0))
        return p₁ * p₂
    else
        p₁ = LaurentPolynomial([pr/(pr-pl),T(1)/(pl-pr)], 0, false, T(0))
        p₂ = LaurentPolynomial([pl/(pl-pr),T(1)/(pr-pl)], 0, false, T(0))
        pⱼ₋₁ = Bubble_laurent(mesh, i, j-1, T)
        return pⱼ₋₁*(p₂ - p₁)
    end
end

function Bubble(mesh::OneDMesh, i::Int, j, T::Type = Float64)
    pⱼ = Bubble_laurent(mesh, i, j, T)
    PiecewiseLaurentPolynomial(mesh, [pⱼ], [i], T(0))
end

function BubbleBasis(mesh::OneDMesh{TM}, T::Type = Float64; order::Int = 1, left::Bool = true, right::Bool = left) where TM
    @assert order ≥ 1
    basis = PiecewiseLaurentPolynomial{T,TM}[]
    for i ∈ eachindex(mesh)[begin:end-1]
        if i ≠ firstindex(mesh) || left
            push!(basis, HatFunctionP1(mesh, i, T))
        end
        for p ∈ 2:order
            push!(basis, Bubble(mesh, i, p+1, T))
        end
    end
    if right
        push!(basis, HatFunctionP1(mesh, lastindex(mesh), T))
    end
    LaurentPolynomialBasis(basis)
end


########################################################################################
#                                   Legendre Basis
########################################################################################

# legendre polynomial rescaled on [m[i], m[i+1]]
@memoize function IntLegendre(mesh::OneDMesh, i::Int, n::Int, T::Type = Float64)
    @assert i ≤ lastindex(mesh) - 1
    pₙ = Legendre(mesh[i], mesh[i+1], n, T, true)
    int_pₙ = integrate(pₙ)
    int_pₙ - int_pₙ(mesh[i])
end

function IntLegendre_element(mesh::OneDMesh, i::Int, n::Int, T::Type = Float64)
    pₙ = IntLegendre(mesh, i, n, T)
    PiecewiseLaurentPolynomial(mesh, [pₙ], [i], T(0))
end

function IntLegendreBasis(mesh::OneDMesh{TM}, T::Type = Float64; order::Int = 1, left::Bool = true, right::Bool = left) where TM
    @assert order ≥ 1
    basis = PiecewiseLaurentPolynomial{T,TM}[]
    for i ∈ eachindex(mesh)[begin:end-1]
        if i != firstindex(mesh) || left       
            push!(basis, HatFunctionP1(mesh, i, T))
        end
        for p ∈ 2:order 
            push!(basis, IntLegendre_element(mesh, i, p, T))
        end
    end
    if right
        push!(basis, HatFunctionP1(mesh, lastindex(mesh), T))
    end
    LaurentPolynomialBasis(basis)
end
