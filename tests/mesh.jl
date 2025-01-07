@inline function findindex(m, x)
    for i in eachindex(m)
        if m[i] > x
            return i-1
        end
    end
    if x == m[end]
        return lastindex(m)
    else
        return lastindex(m)+1
    end
end

@inline function findindex2(m, x)
    idx = searchsortedfirst(m, x)
    return idx > 1 && m[idx - 1] == x ? idx - 1 : idx
end

m = LinRange(0,10,100)

@show findindex(m, 0.1)
@show findindex2(m, 0.1)