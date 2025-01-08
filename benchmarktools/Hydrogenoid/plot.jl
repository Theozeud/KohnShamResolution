# Convergence Plot with Nmesh

function convergence_plot_Nmesh(sols::HydrogenoidConvergenceNmesh, nums = first(sols.num))
    plt = plot( size = (1300,1000), margin = 0.5Plots.cm, legend = :bottomleft , xaxis=:log, yaxis=:log,
                legendfontsize  = 12,  
                titlefontsize   = 12,
                guidefontsize   = 12,
                tickfontsize    = 12)
    xlabel!(plt, "Nmesh")
    ylabel!(plt, "Error")
    title!(plt, "Convergence Plot")
    for prob ∈ sols.probs
        for num ∈ nums
            pos = findfirst(x->x==num, nums)
            plot!(plt, sols.vecNmesh, sols.ΔΛ[prob.name][:,pos], lw = 4, label = prob.name*"-$num", markershape = :x, markersize = 10)
        end
    end
    plt
end 

# Convergence Plot with Rmax

function convergence_plot_Rmax(sols::HydrogenoidConvergenceNmesh, nums = first(sols.num))
    plt = plot( size = (1300,1000), margin = 0.5Plots.cm, legend = legend ? :bottomleft : false, xaxis=:log, yaxis=:log,
                legendfontsize  = 12,  
                titlefontsize   = 12,
                guidefontsize   = 12,
                tickfontsize    = 12)
    xlabel!(plt, "Rmax")
    ylabel!(plt, "Error")
    title!(plt, "Convergence Plot")
    for prob ∈ sols.probs
        for num ∈ nums
            plot!(plt, sols.vecRmax, sols.ΔΛ[prob.name][:,num], lw = 4, label = prob.name*"-$num", markershape = :x, markersize = 10)
        end
    end
    plt
end 