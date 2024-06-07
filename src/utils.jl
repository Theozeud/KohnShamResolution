struct NothingFunction <: Function end
(::NothingFunction)(args...;kwargs...) = nothing

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end