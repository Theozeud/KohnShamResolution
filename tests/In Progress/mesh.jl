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
    if x ≤ m[end]
        return searchsortedlast(m, x)
    else
        return lastindex(m)+1
    end
end

m = LinRange(0,10,600)

@show @btime findindex(m, 5)
@show @btime findindex2(m, 5)