using Plots

# Définir une fonction de potentiel avec une variation rapide
function potential(x)
    return x *exp(-(x - 0.5)^2)  # Puits de potentiel centré à x = 0.5
end

# Créer une grille adaptative
function adaptive_grid(mesh, n, f)
    diff = rand(length(mesh)-1)
    for i in eachindex(diff)
        diff[i] = abs(f(mesh[i+1]) - f(mesh[i]))
    end
    while length(mesh) < n
        argm = argmax(diff)
        mid = (mesh[argm] + mesh[argm+1])/2
        insert!(mesh, argm+1, mid)
        deleteat!(diff, argm)
        insert!(diff, argm, abs(f(mesh[argm+1]) - f(mesh[argm])))
        insert!(diff, argm+1, abs(f(mesh[argm+2]) - f(mesh[argm+1])))
    end
    return mesh
end

function random_adaptive_grid(mesh, n, f)
    diff = rand(length(mesh)-1)
    for i in eachindex(diff)
        diff[i] = abs(f(mesh[i+1]) - f(mesh[i]))
    end
    while length(mesh) < n
        argm1 = rand(eachindex(diff))
        argm2 = rand(eachindex(diff))
        mid1 = (mesh[argm1] + mesh[argm1+1])/2
        mid2 = (mesh[argm2] + mesh[argm2+1])/2
        
        (argm, mid) = if abs(f(mid1) - f(mesh[argm1]))  + abs(f(mesh[argm1+1]) - f(mid1)) >
            abs(f(mid2) - f(mesh[argm2]))  + abs(f(mesh[argm2+1]) - f(mid2))
            (argm1, mid1)
        else
            (argm2, mid2)
        end
        insert!(mesh, argm+1, mid)
        deleteat!(diff, argm)
        insert!(diff, argm, abs(f(mesh[argm+1]) - f(mesh[argm])))
        insert!(diff, argm+1, abs(f(mesh[argm+2]) - f(mesh[argm+1])))
    end
    return mesh
end

function adaptive_grid2(mesh, n, f)
    diff = rand(length(mesh)-1)
    for i in eachindex(diff)
        diff[i] = abs(f(mesh[i+1]) - f(mesh[i]))
    end
    while length(mesh) < n
        mid = (mesh[1] + mesh[2])/2
        argm = 1
        v  = abs(f(mid) - f(mesh[1]))  + abs(f(mesh[2]) - f(mid))
        for i ∈ eachindex(diff)[2:end-1]
            _mid = (mesh[i] + mesh[i+1])/2
            _v = abs(f(_mid) - f(mesh[i]))  + abs(f(mesh[i+1]) - f(_mid))
            if _v > v
                argm = i
                v = _v
                mid = _mid
            end
        end 
        insert!(mesh, argm+1, mid)
        deleteat!(diff, argm)
        insert!(diff, argm, abs(f(mesh[argm+1]) - f(mesh[argm])))
        insert!(diff, argm+1, abs(f(mesh[argm+2]) - f(mesh[argm+1])))
    end

    return mesh
end



mesh = collect(LinRange(0,5,11))
n = 40

x = LinRange(0,5,100)
plot(x, potential.(x), label="true", ls = :dash, lc = :black, lw = 3)

x_adaptive = adaptive_grid2(mesh, n, potential)
plot!(x_adaptive, potential.(x_adaptive), label="grille adaptative", lw = 3)

x_uniforme = LinRange(0,5,15)
plot!(x_uniforme, potential.(x_uniforme), label="grille uniforme", lw = 3)

