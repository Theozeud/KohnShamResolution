function plot_stopping_criteria(sols, tol; title, label)
    # Plot Criteria on Density
    plt_criteria = plot(size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis=:log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_criteria, "Iteration")
    ylabel!(plt_criteria, "Convergence criteria")
    title!(plt_criteria, title)
    lmax = 0
    for (i,sol) ∈ enumerate(sols)
        plot!(plt_criteria, sol.crit, lw = 4, label = label[i], markershape = :x, markersize = 10)
        if length(sol.crit)>lmax
            lmax = length(sol.crit)
        end
    end
    plot!(ones(lmax)*tol, label = "tolerance", ls = :dash, color = :black, lw =4 )

end

#=
    # Plot Criteria on Energy
    plt_energy = plot(  size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis = :log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_energy, "Iteration")
    ylabel!(plt_energy, "Convergence of the Energy")
    title!(plt_energy, title)
    for (i,sol) ∈ enumerate(sols)
        plot!(plt_energy, eachindex(sol.crit)[2:end], 1e-16 .+ abs.(sol.Energyhisto[2:end] .- sol.Energyhisto[1:end-1]), lw = 3, label = label[i], markershape = :x, markersize = 10)
    end
    plot(plt_criteria, plt_energy, layout = (1,2), size = (1100,400))
    =#