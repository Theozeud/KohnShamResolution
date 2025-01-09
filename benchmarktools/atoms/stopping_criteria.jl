function plot_stopping_criteria(sols)
    plt = plot( size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis=:log,
                legendfontsize  = 12,  
                titlefontsize   = 12,
                guidefontsize   = 12,
                tickfontsize    = 12)
    xlabel!(plt, "Iteration")
    ylabel!(plt, "Scf Tolerance")
    title!(plt, "Scf Tolerance through iterations")
    lmax = 0
    for i ∈ eachindex(sols)
        @unpack stopping_criteria_log = sols[i].log
        plot!(plt, stopping_criteria_log, lw = 4, label = sols[i].name, markershape = :x, markersize = 10)
        if length(stopping_criteria_log) > lmax
            lmax = length(stopping_criteria_log)
        end
    end
    #plot!(ones(lmax)*sols[1]., label = "tolerance", ls = :dash, color = :black, lw =4 )
    plt
end