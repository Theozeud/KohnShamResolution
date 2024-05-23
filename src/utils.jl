struct NothingFunction <: Function end
(::NothingFunction)(args...;kwargs...) = nothing