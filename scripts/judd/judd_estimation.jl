using DrWatson
@quickactivate "FractalDimension"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot
using ProgressMeter, Statistics
using JLD2

include("judd.jl")

sample_size = 100
realizations = 100
deg_a = 1
F = 8.

datasets[3] = :lorenz96_chaotic
labels[3] = "Lorenz96 (\$D\$=8)"

datasets[1] = :henon_chaotic
labels[1] = "Hénon map"

datasets[2] = :lorenz_chaotic
labels[2] = "Lorenz system"

function run_exs(sample_size, realizations, F, deg_a)

    name = savename("judd_sensitivity", @dict(sample_size, realizations,deg_a,F),"jld2")

    isfile(datadir("judd", name)) && return

    systems = ["Henon", "Lorenz", "L96-D=5","L96-D=6","L96-D=7","L96-D=8","L96-D=9"]
    systs = [Systems.henon(), Systems.lorenz(),Systems.lorenz96(5, F=F),Systems.lorenz96(6, F=F),Systems.lorenz96(7, F=F),Systems.lorenz96(8, F=F),Systems.lorenz96(9, F=F)]

    Deltas_ic = []

    println("Estimation of sensitivity on initial conditions with fixed sample")
    for (system, ds) in zip(systems, systs)

        println("System: $system")

        ts = regularize(trajectory(ds,Ttr = 100, 5000))
        x = interpoint_distances(ts, sample_size)
        εs, bs, ε₀ = bin_judd(x)

        Deltas_initial_cond = zeros(realizations)
        @showprogress for k in 1:realizations
            d₀, a₀ = rand()*dimension(ds)*1.5, rand(deg_a+1)*500
            d, a = nested_optimizer(εs, bs; ε₀, d₀, a₀, η = 1e-4)
            Deltas_initial_cond[k] = d
        end

        println("For 100 different random initial conditions, the fractal dimension of $system is estimated as $(mean(Deltas_initial_cond))±$(std(Deltas_initial_cond))")
        push!(Deltas_ic, Deltas_initial_cond)
    end


    Deltas_s = []

    println("Estimation of sensitivity on sample with fixed initial conditions")

    for (system, ds) in zip(systems, systs)

        ts = regularize(trajectory(ds,Ttr = 100, 5000))

        d₀, a₀ = rand()*dimension(ds)*1.5, rand(deg_a+1)*500
        Deltas_samples = zeros(realizations)
        @showprogress for k in 1:realizations
            x = interpoint_distances(ts, sample_size)
            εs, bs, ε₀ = bin_judd(x)
            d, a = nested_optimizer(εs, bs; ε₀, d₀, a₀, η = 1e-4)
            Deltas_samples[k] = d
        end

        println("For 100 different random samples, the fractal dimension of $system is estimated as $(mean(Deltas_samples))±$(std(Deltas_samples))")
        push!(Deltas_s, Deltas_samples)
    end

    result = @strdict systems Deltas_ic Deltas_s
    @tagsave(datadir("judd", name), result)
end

function plot_exs(sample_size,realizations,F,deg_a)
    name = savename("judd_sensitivity", @dict(sample_size, realizations,deg_a,F),"jld2")
    plotname = savename("judd_sensitivity", @dict(sample_size, realizations,deg_a,F),"png")

    # isfile(plotsdir("judd", plotname)) && return

    @load datadir("judd", name) systems Deltas_ic Deltas_s

    popat!(Deltas_s,findall(x->x=="Lorenz",systems)[1])
    popat!(systems,findall(x->x=="Lorenz",systems)[1])

    fig = figure()

    pos = collect(1:length(Deltas_s))
    fp = Dict("marker"=>"o", "markerfacecolor"=>COLORS[1], "markersize"=>6,
                  "linestyle"=>LINESTYLES[1], "markeredgecolor"=>COLORS[1])
    box1 = boxplot(Deltas_s, notch = true, positions = pos, flierprops=fp,patch_artist=true)
    for item in ["boxes", "whiskers", "fliers", "medians", "caps"]
        setp(box1[item], color=COLORS[1])
    end

    xticks(pos, systems, rotation=45)
    xlabel("System")
    ylabel(L"$ \Delta^{(J)}$")
    tight_layout()
    wsave(plotsdir("judd", plotname), fig)
end

run_exs(sample_size,realizations,F,deg_a)

plot_exs(sample_size,realizations,F,deg_a)
