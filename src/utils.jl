struct NothingFunction <: Function end
(::NothingFunction)(args...;kwargs...) = nothing

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

function intersection_with_indices(A, B)
    # Find intersection of A and B
    intersect_elements = intersect(A, B)

    # Create a vector of tuples (element, index_in_A, index_in_B)
    return  [(indfirst(x -> x == el, A), findfirst(x -> x == el, B)) for el in intersect_elements]
end