using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot
using ProgressMeter, Statistics
using JLD2

include("koch.jl")

function run_exs_takens(k, F, N, eps_step)

    name = savename("takens_epsilonmax_regularized", @dict(k, N,F,eps_step),"jld2")
    isfile(datadir("takens", name)) && return

    systems = ["Koch", "Henon", "Lorenz", "L96_D5","L96_D6","L96_D7","L96_D8","L96_D9"]
    systs = [pointskoch(points, k), Systems.henon(), Systems.lorenz(),Systems.lorenz96(5, F=F),Systems.lorenz96(6, F=F),Systems.lorenz96(7, F=F),Systems.lorenz96(8, F=F),Systems.lorenz96(9, F=F)]

    Delta_T = Dict()

    for (system, ds) in zip(systems, systs)

        println("System: $system")

        X = nothing

        if system == "Koch"
            points = [[0.0; 0.0], [1.0; 0.0]]
            koch = pointskoch(points, k)
            # making sure we have the same amount of points for the Koch curve
            pt = rand(1:length(koch), 40000)
            X = regularize(Dataset(koch[pt]))
        else
            dt = system == "Henon" ? 1 : 0.1
            T = N * dt
            ts = trajectory(ds, T; dt = dt, Ttr = 100)
            X = regularize(ts)
        end

        εs = 0.2:eps_step:1
        α, αmin, αmax = zeros(length(εs)), zeros(length(εs)), zeros(length(εs))

        @showprogress for (i,ε) in enumerate(εs)
            try
                α[i], αmax[i], αmin[i] = takens_best_estimate(X, ε)
            catch
                # An error above means that the estimator diverged
                α[i], αmax[i], αmin[i] = NaN, NaN, NaN
            end
        end
        result = @strdict α αmin αmax

        Delta_T[system] = result
    end
    @tagsave(datadir("takens", name), Delta_T)
end

function plot_ex_takens(F, N, eps_step,k)

    name = savename("takens_epsilonmax_regularized", @dict(N,F,eps_step,k),"jld2")
    plotname = savename("takens_epsilonmax_regularized_all", @dict(N,F, eps_step,k),"png")

    systems = ["Koch", "Henon", "L96_D5","L96_D7","L96_D9"]
    labels = Dict("Koch"=> "Koch", "Henon"=> "Hénon", "L96_D5"=>"Lorenz-96, \$D = 5\$","L96_D7"=>"Lorenz-96, \$D = 7\$","L96_D9"=>"Lorenz-96, \$D = 9\$")

    @load datadir("takens", name) Koch Henon L96_D5 L96_D7 L96_D9

    system_results = [Koch, Henon, L96_D5, L96_D7, L96_D9]

    @load datadir("lorenz96_Dmax=13_F=$(Int(F))_T=100000_ics=10.jld2") K

    K2 = [log(4)/log(3), 1.26] # analytial solution for Koch, Δ^(L) for Hénon

    Ks = [K2...,mean.(K[[1,3,5]])...]

    εs = 0.2:eps_step:1

    fig, axs = plt.subplots(2, sharex=true,  gridspec_kw=Dict("hspace"=> .1),constrained_layout=false)

    mylinestyles = [":", "--", "-.",":", "--"]

    for (k,(system,system_result)) in enumerate(zip(systems, system_results))

        ax = (system == "Henon" || system == "Koch") ? axs[2] : axs[1]

        ax.fill_between(εs, system_result["αmin"], system_result["αmax"], color = COLORS[k], alpha = .5)
        ax.plot(εs,system_result["α"], color = COLORS[k],linestyle = mylinestyles[k],label = labels[system])

        ax.plot(εs, ones(length(εs))*Ks[k], color = COLORS[k], linestyle = "-", lw = 3)
        ax.set_xlim(.2,1)
    end

    handles, labels = axs[1].get_legend_handles_labels()
    handles2, labels2 = axs[2].get_legend_handles_labels()
    xlabel(L"$\varepsilon_{\mathrm{max}}$")
    fig.text(0.04, 0.55, L"$\Delta^{(T)}$", va="center", rotation="vertical", fontsize=30)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.13, top = .84)
    fig.legend([handles...,handles2...],[labels...,labels2...], loc = "upper center", ncol = 3)

    tight_layout()
    wsave(plotsdir("takens", plotname), fig)
    close("all")
end


N = 40000
F = 8.
k = 8
run_exs_takens(k,F,N,.01)
plot_ex_takens(F, N,.01,k)
