
function convergenceNmesh(vecNmesh::AbstractVector, problems; nums = [1])
    Errors = Dict()
    for problem ∈ problems
        @unpack T, name = problem
        println(name)
        ϵ = zeros(T, length(vecNmesh), length(nums))
        @inbounds for i ∈ eachindex(vecNmesh)
            newprob = AtomProblem(problem; Nmesh = vecNmesh[i])
            @time "Nmesh = $(vecNmesh[i])" sol = groundstate(newprob)
            ϵ[i,:] .= sol.orbitals_energy[nums]
        end
        error = zeros(T, length(vecNmesh)-1, length(nums))
        for i ∈ axes(error,1)
            error[i,:] = abs.(ϵ[i,:] .- ϵ[end,:])
        end
        Errors[name] = error
    end
    AtomConvergenceNmesh(vecNmesh, Errors, nums)
end


function convergencePlotNmesh(sols::AtomConvergenceNmesh, nums = first(sols.num))
    plt = plot( size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis=:log,
                legendfontsize  = 12,  
                titlefontsize   = 12,
                guidefontsize   = 12,
                tickfontsize    = 12)
    xlabel!(plt, "Nmesh")
    ylabel!(plt, "Error")
    title!(plt, "Convergence Plot")
    for key ∈ keys(sols.Error)
        for num ∈ nums
            plot!(plt, sols.vecNmesh[1:end-1], sols.Error[key][:,num], lw = 4, label = key*"-$num", markershape = :x, markersize = 10)
        end
    end
    plt
end

function convergenceRmax(vecRmax::AbstractVector, problems; nums = [1])
    Errors = Dict()
    for problem ∈ problems
        @unpack T, name = problem
        println(name)
        ϵ = zeros(T, length(vecRmax), length(nums))
        @inbounds for i ∈ eachindex(vecRmax)
            newprob = AtomProblem(problem; Rmax = vecRmax[i])
            @time "Rmax = $(vecRmax[i])" sol = groundstate(newprob)
            ϵ[i,:] .= sol.orbitals_energy[nums]
        end
        error = zeros(T, length(vecRmax)-1, length(nums))
        for i ∈ axes(error,1)
            error[i,:] = abs.(ϵ[i,:] .- ϵ[end,:])
        end
        Errors[name] = error
    end
    AtomConvergenceRmax(vecRmax, Errors, nums)
end


function convergencePlotRmax(sols::AtomConvergenceRmax, nums = first(sols.num))
    plt = plot( size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis=:log,
                legendfontsize  = 12,  
                titlefontsize   = 12,
                guidefontsize   = 12,
                tickfontsize    = 12)
    xlabel!(plt, "Rmax")
    ylabel!(plt, "Error")
    title!(plt, "Convergence Plot")
    for key ∈ keys(sols.Error)
        for num ∈ nums
            plot!(plt, sols.vecRmax[1:end-1], sols.Error[key][:,num], lw = 4, label = key*"-$num", markershape = :x, markersize = 10)
        end
    end
    plt
end
