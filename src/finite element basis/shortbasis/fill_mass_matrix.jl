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
    A[1,1] = p1.left ? (T(mesh[I+1]) - T(mesh[I]))/T(3) : T(2) * (T(mesh[I+1]) - T(mesh[I]))/T(3)
    for I ∈ axes(A,1)[2:end]
        A[I, I] = T(2) * (T(mesh[I+1]) - T(mesh[I]))/T(3)
        A[I, I-1] = (T(mesh[I+1]) - T(mesh[I]))/T(6)
        A[I-1, I] = A[I, I-1]
    end
    if p1.right
        A[end, end] = (T(mesh[end]) - T(mesh[end-1]))/T(3)
    end
    nothing
end


########################################################################################
#                    Integrated Legendre Polynomials mass matrix
########################################################################################
function fill_mass_matrix!(ilb::IntLegendreElements{true}, ::OneDMesh, A)
    T = eltype(ilb)
    for I ∈ axis(A,1)
        A[I,I] = T(1)
    end
    nothing
end

function fill_mass_matrix!(ilb::IntLegendreElements{false}, mesh::OneDMesh, A)
    for I ∈ axis(A,1)
        A[I, I] = (T(mesh[I+1]) - T(mesh[I]))/(T(ilb.bsup) - T(ilb.binf))
    end
    nothing
end

