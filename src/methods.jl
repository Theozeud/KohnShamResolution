struct Method1 <: AbstractKohnShamResolutionMethod end

for_sphericalsymmetry(::Method1) = true
for_cylindricalsymmetry(::Method1) = false


struct CacheMethod1 <: AbstractKohnShamCache
    A
    M₀
    M₋₁
    M₋₂
end

function init_cache(::Method1)

    # Init matrix

end


function performstep!(method, solver)
    # step 1 : find potential 
    
    # step 2 : compute an approximation of the exchange correlation term

    # step 3 : résolution du problème aux valeurs propres blocs par blocs

    # step 4 : building a certain matrix

    # step 5 : update this matrix with a convex approach

end

function stopping_criteria(method, solver)


end