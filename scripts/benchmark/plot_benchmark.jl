using DrWatson
@quickactivate "FractalDimension"
include(srcdir("style.jl"))
using PyPlot, DataFrames, CSVFiles


# %% computation times for different lengths

data = wload(datadir("benchmark", "benchmark_number_NofN=10_base=10_lowerN=3_upperN=5.csv")) |> DataFrame
fig, axs = subplots(2,1; sharex = true, figsize = (10, 10))
for (i, name) in enumerate(unique(data.Method))
    color = "C$(i-1)"

    lorenz_data = data[(data.Model .== "lorenz") .& (data.Method .== name), :]
    axs[1].plot(lorenz_data.N, lorenz_data.Time, label = name, color = color)

    hénon_data = data[(data.Model .== "henon") .& (data.Method .== name), :]
    axs[2].plot(hénon_data.N, hénon_data.Time, color = color)
end
axs[2].set_xlabel("\$N\$")
axs[2].set_ylabel("\$t_{Hénon}[\\mathrm{ns}]\$")
axs[1].set_ylabel("\$t_{Lorenz63}[\\mathrm{ns}]\$")

axs[1].set_yscale("log")
axs[2].set_yscale("log")
axs[1].set_xscale("log")
axs[2].set_xscale("log")
for ax in axs; ax.legend(loc = "upper left"); ax.grid(); end
axs[1].tick_params(labelbottom=false)
axs[1].set_title("minimum computation time")
wsave(plotsdir("benchmarks", "comparing_btime_N.png"), fig)

# %% computation times for different system dimensions

data = wload(datadir("benchmark", "benchmark_dim_N=10000_ddim=1_lowerdim=4_upperdim=16.csv")) |> DataFrame
fig, axs = subplots(1,1; sharex = true, figsize = (10, 16))
for (i, method) in enumerate(unique(data.Method))
    color = "C$(i-1)"

    x = data[data.Method .== method, :Dimension]
    y = data[data.Method .== method, :Time]
    axs.plot(x, y, color = color, label = method)
end
axs.set_xlabel("\$\\mathrm{Dimension}\$")
axs.set_ylabel("\$t [\\mathrm{ns}]\$")

axs.set_yscale("log")

axs.legend(loc = "upper left")
axs.grid()

axs.set_title("minimum computation time (N=10000)")
wsave(plotsdir("benchmarks", "comparing_btime_dim.png"), fig)

# %% comparison of computation times for length and dimensionality

path1 = "benchmark_dim_N=50000_ddim=1_lowerdim=4_upperdim=16.csv"
path2 = "benchmark_number_NofN=10_base=10_lowerN=3_upperN=5.csv"
titlename = "\$ \\mathrm{Computation}\\;\\mathrm{Time}\\;\\mathrm{for}\\;\\mathrm{Calculation}\\;\\mathrm{of}\\;\\mathrm{Probabilities} \$"
labels = [
    "\$ C_2^\\mathrm{classic} \$",
    "\$ C_2^\\mathrm{tree} \$",
    "\$ C_2^\\mathrm{box} \$",
    "\$ C_2^\\mathrm{prism} \$",
    "\$ H_1^\\mathrm{Molteno} \$",
    "\$ H_1^\\mathrm{classic} \$",
]

function plot_dim_length(title, labels, path1, path2; save_figure = true)
    data_dim = wload(
        datadir(
                "benchmark",
                path1,
            )) |> DataFrame
    data_var = wload(
            datadir(
                "benchmark",
                path2,
            )) |> DataFrame
    data_length = data_var[(data_var.Model .== "lorenz"), :]

    fig, axs = subplots(2,1; sharex = false)

    c_names = Dict()
    ls_names = Dict()
    agents = []
    for (i, method) in enumerate(unique(data_dim.Method))
        color = "C$(i-1)"
        ls = LINESTYLES[i]
        c_names[method] = color
        ls_names[method] = ls

        x = data_dim[data_dim.Method .== method, :Dimension]
        y = data_dim[data_dim.Method .== method, :Time]
        agent, = axs[1].plot(x, y, color = color, label = method, ls = ls)
        push!(agents, agent)
    end
    for method in unique(data_length.Method)
        if method == "Correlation"
            color = c_names["Classic Correlation"]
            ls = ls_names["Classic Correlation"]
        elseif method == "Generalized Dimension"
            color = c_names["Generalized Entropy"]
            ls = ls_names["Generalized Entropy"]
        else
            color = c_names[method]
            ls = ls_names[method]
        end

        x = data_length[data_length.Method .== method, :N]
        y = data_length[data_length.Method .== method, :Time]
        axs[2].plot(x, y, color = color, label = nothing, ls = ls)
    end

    for ax in axs; ax.set_yscale("log"); ax.grid(true); end

    axs[1].set_xlabel("\$ \\mathrm{Dimension} \$")
    axs[1].set_xticks(4:4:16)
    axs[1].set_xlim((4,16))
    axs[1].set_title("\$ \\mathrm{Lorenz-96} \\quad N=50000 \$")

    axs[1].set_ylabel("\$t [\\mathrm{ns}]\$")
    axs[1].set_yticks(10 .^ (7:12))
    leg = axs[1].legend(
        handles = agents,
        labels = labels,
        bbox_to_anchor=(0., 1.25, 1., .102), loc="lower left",
        ncol=3, mode="expand", borderaxespad=0, handlelength=2,
        title = title,
    )


    axs[2].set_xlabel("\$\\mathrm{N}\$")
    axs[2].set_yticks(10 .^ (5:2:12))
    axs[2].set_xlim((10^3, 10^5))
    axs[2].set_xscale("log")
    axs[2].set_title("\$ \\mathrm{Lorenz-63} \$")
    #axs[2].legend().remove()

    axs[2].set_ylabel("\$t [\\mathrm{ns}]\$")

    axs[1].add_artist(leg)
    fig.tight_layout()
    fig.subplots_adjust(top = 1.0, hspace = 0.5, left = 0.2)

    fig
end

fig = plot_dim_length(title, labels, path1, path2)
fig.tight_layout(;pad = 0.3)
wsave(plotsdir("benchmarks", "comparing_btime_dim_N.png"), fig)
