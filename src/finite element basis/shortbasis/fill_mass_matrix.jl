########################################################################################
#                                   P1 mass matrix
########################################################################################
function fill_mass_matrix!(p1::P1Elements{true}, ::OneDMesh, A)
    if p1.bsup != 1 || p1.binf != -1
        @warn "Mass matrix is filled like if the bounds of P1 elements are scaled on [-1,1]!"
    end
    T = eltype(p1)
    A[1,1] = p1.left ? T(0.5) : T(1)
    for I ∈ axes(A,1)[2:end]
        A[I,I] = T(1)
        A[I, I-1] = T(0.25)
        A[I-1, I] = A[I, I-1]
    end
    if p1.right
        A[end, end] = T(0.5)
    end
    nothing
end

function fill_mass_matrix!(p1::P1Elements{false}, mesh::OneDMesh, A)
    if p1.bsup != 1 || p1.binf != -1
        @warn "Mass matrix is filled like if the bounds of P1 elements are scaled on [-1,1]!"
    end
    T = eltype(p1)
    A[1,1] = p1.left ? (T(mesh[2]) - T(mesh[1]))/T(3) : (T(mesh[3]) - T(mesh[1]))/T(3)
    for I ∈ axes(A,1)[2:end-1]
        A[I, I] = (T(mesh[I+1+ !p1.left]) - T(mesh[I-1+ !p1.left]))/T(3)
        A[I, I-1] = (T(mesh[I+ !p1.left]) - T(mesh[I-1+ !p1.left]))/T(6)
        A[I-1, I] = A[I, I-1]
    end
    A[end, end-1] = (T(mesh[end-1+p1.right]) - T(mesh[end-2+p1.right]))/T(6)
    A[end-1, end] = A[end, end-1]
    A[end, end] = p1.right ? (T(mesh[end]) - T(mesh[end-1]))/T(3) : (T(mesh[end]) - T(mesh[end-2]))/T(3)

    nothing
end


########################################################################################
#                    Integrated Legendre Polynomials mass matrix
########################################################################################
function fill_mass_matrix!(ilb::IntLegendreElements{true}, ::OneDMesh, A)
    T = eltype(ilb)
    for I ∈ axes(A,1)
        A[I,I] = T(1)
    end
    nothing
end

function fill_mass_matrix!(ilb::IntLegendreElements{false}, mesh::OneDMesh, A)
    T = eltype(ilb)
    for i ∈ eachindex(mesh)[1:end-1]
        for n ∈ 1:ilb.size
            I = (i-1)*ilb.size + n
            A[I,I] = (T(mesh[i+1]) - T(mesh[i]))/(T(ilb.bsup) - T(ilb.binf))
        end
    end
    nothing
end


########################################################################################
#                    Difference Legendre Polynomials mass matrix
########################################################################################
function fill_mass_matrix!(dlb::DiffLegendreElements{true}, ::OneDMesh, A)
    T = eltype(dlb)
    C(j) = 1/sqrt(2) * sqrt((2j+3)*(2j-1)/(4j+2))
    for I ∈ axes(A,1)
        A[I,I] = T(1)
        if I ≥ 2
            A[I-2, I] = A[I, I-2]
        end
        if I ≤ axes(A,1) - 2
            A[I, I+2] = -T(2)/T(2*I+3) * C(I) * C(I+2)
        end
    end
    nothing
end