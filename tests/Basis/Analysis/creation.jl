using KohnShamResolution

# Creation of elements
T = Float64
ordermin = 2
ordermax = 3

Nmesh = 1000
Rmax  = 20
mesh = linmesh(0, Rmax, Nmesh; T = T)

# Tests integrals on elements
function analysis_ShortIntLegendreBasis(mesh, T; Rcut::Real = last(mesh), kwargs...)
    println("Creation elements")
    @time intlegelement = IntLegendreElements(T; kwargs...)
    Ncut = min(KohnShamResolution.findindex(mesh, Rcut), lastindex(mesh))
    size = intlegelement.size * (Ncut - 1)
    println("Initialisation InfoElements")
    @time infos = Vector{InfoElement{T}}(undef, size)
    println("Fill InfoElements")
    @time for i ∈ eachindex(mesh)[1:Ncut-1]
        segments = [i]
        shifts = [KohnShamResolution.shift(T, mesh[i], mesh[i+1], intlegelement.binf, intlegelement.bsup)]
        invshifts = [KohnShamResolution.shift(T, intlegelement.binf, intlegelement.bsup, mesh[i], mesh[i+1])]
        for n ∈ 1:intlegelement.size
            index = [n]
            infos[(i-1) * intlegelement.size + n] = InfoElement(index, segments, shifts, invshifts)
        end
    end
    println("Creation of the basis")
    analysis_ShortPolynomialBasis(intlegelement, mesh, size, infos)
end


function analysis_ShortPolynomialBasis(elements, mesh::Mesh, size::Int, infos) 
        
    # Basics checks for consistency
    @assert size == length(infos)
    for info ∈ infos
        @assert KohnShamResolution.getsegments(info) ⊆ eachindex(mesh)[1:end-1]
    end

    # Creating the coupling indices for matrix
    println("Creation of coupling_index 2")
    @time begin
    coupling_index = CartesianIndex{2}[]
    @inbounds for i in eachindex(infos)
        @inbounds for j in eachindex(infos)[i:end]
            if !isdisjoint(KohnShamResolution.getsegments(infos[i]), KohnShamResolution.getsegments(infos[j]))
                push!(coupling_index, CartesianIndex(i, j))
            else
                break
            end
        end
    end end

    # Creating the coupling indices for 3-tensor
    println("Creation of coupling_index 3")
    @time begin
    coupling_index3 = CartesianIndex{3}[]
    @inbounds for I ∈ coupling_index
        i = I[1]
        j = I[2]
        S = intersect(KohnShamResolution.getsegments(infos[i]), KohnShamResolution.getsegments(infos[j]))
        @inbounds for k ∈ eachindex(infos)[i:end]
            if !isdisjoint(S, KohnShamResolution.getsegments(infos[k]))
                push!(coupling_index3, CartesianIndex(i, j, k))
            else
                break
            end
        end
    end end

end


analysis_ShortIntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax)
nothing


### Combined Basis

function analysis_ShortP1IntLegendreBasis(mesh::Mesh, T::Type = Float64; left = false, right = false, ordermin = 2, ordermax = 2, first = true, Rcut = last(mesh))
    if ordermax == 1
        return ShortP1Basis(mesh, T; left = left, right = right)
    else
        intleg = ShortIntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax, Rcut = Rcut)
        p1 = ShortP1Basis(mesh, T; left = left, right = right)
        if first
            return @time analysis_CombineShortPolynomialBasis(p1, intleg)
        else
            return @time analysis_CombineShortPolynomialBasis(intleg, p1)
        end
    end
end  



function analysis_CombineShortPolynomialBasis(basisVector...)
    size = sum([length(basis) for basis ∈ basisVector])
    blocks  = Vector{InfoBlock}(undef, length(basisVector)*(length(basisVector)+1)÷2)
    blocks3 = Vector{InfoBlock}(undef, length(basisVector)*(length(basisVector)+1)*(length(basisVector)+2)÷6)
    size_i = 1
    ib = 1
    ib3 = 1
    cumul_index = zeros(Int, length(basisVector))
    for i ∈ eachindex(basisVector)      # small loop
        cumul_index[i] = size_i
        axes1 = size_i:size_i+length(basisVector[i])-1
        size_j = 1
        for j ∈ 1:i                     # small loop
            axes2 = size_j:size_j+length(basisVector[j])-1
            interaction_index = CartesianIndex{2}[]
            size_k = 1
            for k ∈ 1:j                 # small loop
                axes3 = size_k:size_k+length(basisVector[k])-1
                interaction_index3 = CartesianIndex{3}[]
                for n in eachindex(basisVector[i].infos)
                    for m in eachindex(basisVector[j].infos)
                        intersect_nm = intersect(KohnShamResolution.getsegments(basisVector[i], n), KohnShamResolution.getsegments(basisVector[j], m))
                        if !isempty(intersect_nm )
                            if k == 1
                                push!(interaction_index, CartesianIndex(n, m))
                            end
                            for p in eachindex(basisVector[k].infos)
                                if !isempty(intersect(intersect_nm , KohnShamResolution.getsegments(basisVector[k], p)))
                                    push!(interaction_index3, CartesianIndex(n,m,p))
                                end                                    
                            end
                        else
                            break
                        end
                    end
                end
                size_k += length(basisVector[k])
                block3 = InfoBlock(CartesianIndex(i,j,k), (axes1, axes2, axes3), i == j == k, interaction_index3)
                blocks3[ib3] = block3
                ib3+=1
            end
            size_j += length(basisVector[j])
            block = InfoBlock(CartesianIndex(i,j), (axes1, axes2), i == j, interaction_index)
            blocks[ib] = block
            ib+=1
        end
        size_i += length(basisVector[i])
    end
    blocks, blocks3
end


@time blocks, blocks3 = analysis_ShortP1IntLegendreBasis(mesh, T; ordermax = ordermax)

@time basiss = ShortP1IntLegendreBasis(mesh, T; ordermax = ordermax)

tblock = basis.blocks
tblock3 = basis.blocks3
nothing