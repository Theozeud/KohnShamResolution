using Plots
using PrettyTables

function compute_sol(;z, N, Rmax, Nmesh, lₕ, maxiter, oda, tol = 1e-5, T = Float64)
    method = ConstantODA(T(oda))
    m = linmesh(zero(T), Rmax, Nmesh)
    KM = KohnShamExtended(z = z, N = N)

    basis_p1 = ShortP1Basis(m, T; left = false, right = false, normalize = true)
    basis_intleg2 = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 2, normalize = true)
    basis_intleg3 = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 3, normalize = true)
    basis_intleg4 = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 4, normalize = true)
    basis_intleg5 = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 5, normalize = true)

    # Final Discretization
    D1 = KohnShamRadialDiscretization(lₕ, basis_p1, m)
    D2 = KohnShamRadialDiscretization(lₕ, basis_intleg2, m)
    D3 = KohnShamRadialDiscretization(lₕ, basis_intleg3, m)
    D4 = KohnShamRadialDiscretization(lₕ, basis_intleg4, m)
    D5 = KohnShamRadialDiscretization(lₕ, basis_intleg5, m)

    # Solution
    @time "With P1" sol1 = groundstate(KM, D1, method; tol = tol, hartree = true, maxiter = maxiter)
    @time "With IntLeg2" sol2 = groundstate(KM, D2, method; tol = tol, hartree = true, maxiter = maxiter)
    @time "With IntLeg3" sol3 = groundstate(KM, D3, method; tol = tol, hartree = true, maxiter = maxiter)
    @time "With IntLeg4" sol4 = groundstate(KM, D4, method; tol = tol, hartree = true, maxiter = maxiter)
    @time "With IntLeg5" sol5 = groundstate(KM, D5, method; tol = tol, hartree = true, maxiter = maxiter)
    sols = [sol1, sol2, sol3, sol4, sol5]
    # Title
    title = "Rmax = $Rmax, z = $z, t = $oda, N = $N, Nmesh = $Nmesh"
    # Label
    label = ["P1","IntLeg2","IntLeg3", "IntLeg4", "IntLeg5"]
    # Return
    sols, title, label
end


function plot_crit(sols, tol; title, label)
    plt_criteria = plot( size = (1300,1000), margin = 0.5Plots.cm, legend = :outertopright, yaxis=:log,
    legendfontsize  = 18,  
    titlefontsize   = 18,
    guidefontsize   = 18,
    tickfontsize    = 18)
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
    plt_criteria
end

function plot_Ehisto(sols; title, label)
    plt_Ehisto = plot(  size = (1300,1000), margin = 0.5Plots.cm, legend = :outertopright,
                        legendfontsize  = 11,  
                        titlefontsize   = 11,
                        guidefontsize   = 11,
                        tickfontsize    = 11)
    xlabel!(plt_Ehisto, "Iteration")
    ylabel!(plt_Ehisto, "Evolution of ϵ₁")
    title!(plt_Ehisto, title)
    for (i,sol) ∈ enumerate(sols)
        Ehisto = [min(sol.ϵhisto[i]...) for i∈ eachindex(sol.ϵhisto)]
        plot!(plt_Ehisto, Ehisto, lw = 3, label = label[i], markershape = :x, markersize = 10)
    end
    plt_diffE = plot(  size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis = :log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_diffE, "Iteration")
    ylabel!(plt_diffE, "Convergence of ϵ₁")
    title!(plt_diffE, title)
    for (i,sol) ∈ enumerate(sols)
        Ehisto = [min(sol.ϵhisto[i]...) for i∈ eachindex(sol.ϵhisto)]
        plot!(plt_diffE, eachindex(sol.crit)[2:end], abs.(Ehisto[2:end] .- Ehisto[1:end-1]), lw = 3, label = label[i], markershape = :x, markersize = 10)
    end
    plt_Ehisto, plt_diffE
end

function plot_density(sols, Rmax; title, label)
    plt_ρ = plot(  size = (1300,1000), margin = 0.5Plots.cm, legend = :outertopright,
                        legendfontsize  = 18,  
                        titlefontsize   = 18,
                        guidefontsize   = 18,
                        tickfontsize    = 18)
    xlabel!(plt_ρ, "r")
    ylabel!(plt_ρ, "ρ")
    title!(plt_ρ, title)
    X = range(0, Rmax, 1000)
    for (i,sol) ∈ enumerate(reverse(sols))
        plot!(plt_ρ, X, sol.ρ.(X), lw = 3, label = label[length(sols) - i + 1])
    end
    plt_ρ
end

function plot_Energy(sols; title, label)
    plt = plot(  size = (1300,1000), margin = 1Plots.cm, legend = :outertopright,
                        legendfontsize  = 13,  
                        titlefontsize   = 13,
                        guidefontsize   = 13,
                        tickfontsize    = 13)
    xlabel!(plt, "Iteration")
    ylabel!(plt, "Evolution of the Energy")
    title!(plt, title)
    for (i,sol) ∈ enumerate(sols)
        plot!(plt, sol.Energyhisto, lw = 3, label = label[i], markershape = :x, markersize = 10)
    end

    plt_diff = plot(size = (1300,1000), margin =1Plots.cm, legend = :outertopright, yaxis = :log,
                    legendfontsize  = 18,  
                    titlefontsize   = 18,
                    guidefontsize   = 18,
                    tickfontsize    = 18)
    xlabel!(plt_diff, "Iteration")
    ylabel!(plt_diff, "Convergence of the Energy")
    title!(plt_diff, title)
    for (i,sol) ∈ enumerate(sols)
        plot!(plt_diff, eachindex(sol.crit)[2:end], 1e-16 .+ abs.(sol.Energyhisto[2:end] .- sol.Energyhisto[1:end-1]), lw = 3, label = label[i], markershape = :x, markersize = 10)
    end
    plt_diff
    #plot(plt, plt_diff, layout = (1,2), size = (1650,600))
end



function _groundstate_withNmesh(model, method, order; Rmax, Nmesh, lₕ, maxiter, tol = 1e-5, T = Float64)
    EnergyArray = []
    for N ∈ Nmesh
        sol,_,_ = _groundstate(model, method, order; Rmax = Rmax, Nmesh = N, lₕ = lₕ, maxiter = maxiter, tol = tol,  T = T)
        push!(EnergyArray, sol.Energy)
    end
    plt_diff = plot(size = (1300,1000), margin = 1Plots.cm, legend = :outertopright, yaxis = :log,
                    legendfontsize  = 18,  
                    titlefontsize   = 18,
                    guidefontsize   = 18,
                    tickfontsize    = 18)
    xlabel!(plt_diff, "Nmesh")
    ylabel!(plt_diff, "Convergence of the Energy")
    title!(plt_diff, "Rmax = $Rmax, z = $(model.z), t = $(method.t), N = $(model.N)")
    plot!(plt_diff, Nmesh[2:end], 1e-16 .+ abs.(EnergyArray[2:end] .- EnergyArray[1:end-1]), lw = 3, label = order == 1 ? "P1" : "IntLeg$order", markershape = :x, markersize = 10)
    plt_diff
end

function _groundstate(model, method, order; Rmax, Nmesh, lₕ, maxiter, tol = 1e-5, T = Float64)
    m = linmesh(zero(T), Rmax, Nmesh)
    if order == 1 
        basis = ShortP1Basis(m, T; left = false, right = false, normalize = true)
    else
        basis = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = order, normalize = true)
    end
    # Final Discretization
    D = KohnShamRadialDiscretization(lₕ, basis, m)
    # Solution
    @time sol = groundstate(model, D, method; tol = tol, hartree = true, maxiter = maxiter)

    # Title
    title = "Rmax = $Rmax, z = model.z, oda = method.t, N = model.N, Nmesh = $Nmesh"
    # Label
    label = order == 1 ? "P1" : "IntLeg$order"
    # Return
    sol, title, label
end


function recap_value(sols; title, label)

    



end