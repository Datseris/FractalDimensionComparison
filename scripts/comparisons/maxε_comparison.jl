# Compare the sizes of the estimated `r0` from the Bueno-Orovio or Theiler
# or the standard one from ChaosTools.jl (attractor size)
using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))

datasets = Vector(undef, 6)
labels = Vector{String}(undef, 6)

N = Int(1e5)

datasets[1] = :roessler_chaotic
labels[1] = "Rössler"
datasets[2] = :henonheiles_chaotic
labels[2] = "Hénon-Heiles"
datasets[3] = :lorenz96_chaotic
labels[3] = "Lorenz96 (\$D\$=8)"
datasets[4] = :henon_chaotic
labels[4] = "Hénon map"
datasets[5] = :towel_chaotic
labels[5] = "Towel map"
datasets[6] = :coupled_logistics
labels[6] = "CLM (\$D\$=8)"

fig = figure()
ax = gca()

for i in 1:length(datasets)
    # These are the parameters that change between the different elements of the plot
    data = datasets[i]
    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict N data
    if (data == :lorenz96_chaotic) || (data == :henonheiles_chaotic)
        params["D"] = 8
    end

    data_producing_function = getfield(Data, data)
    X = data_producing_function(; dict2ntuple(params)...)

    P = ChaosTools.autoprismdim(X)
    t0, ε0 = estimate_r0_theiler(X)
    r0, ε0 = estimate_r0_buenoorovio(X, P)

    mi, ma = minmaxima(X)
    d0 = mean(ma - mi)

    x = 2
    a = log(d0/ε0)
    a2 = log(d0/ε0) - x # by default we reduce d0 by x orders of magnitude (base e)
    b = log(r0/ε0)
    c = log(t0/ε0)
    ax.scatter(i, a; marker = "^", color = "C$(i-1)", s = 200, alpha = 0.75, linewidths = 1, edgecolors = "k")
    ax.scatter(i, a2; marker = "o", color = "C$(i-1)", s = 200, alpha = 0.75,linewidths = 1, edgecolors = "k")
    ax.scatter(i, b; marker = "s", color = "C$(i-1)", s = 200, alpha = 0.75,linewidths = 1, edgecolors = "k")
    ax.scatter(i, c; marker = "D", color = "C$(i-1)", s = 200, alpha = 0.75,linewidths = 1, edgecolors = "k")
end

custom_lines = [
    matplotlib.lines.Line2D([0], [0], ls = "None",marker = "^", color = "k", ms = 10),
    matplotlib.lines.Line2D([0], [0], ls = "None",marker = "o", color = "k", ms = 10),
    matplotlib.lines.Line2D([0], [0], ls = "None",marker = "s", color = "k", ms = 10),
    matplotlib.lines.Line2D([0], [0], ls = "None",marker = "D", color = "k", ms = 10),
]
ax.legend(custom_lines, ["mean data length", "m.d.l. / \$e^2\$", "Bueno-Orovio", "Theiler"];
    ncol = 1
)
ax.set_ylabel("\$\\log(\\varepsilon_{max}/\\varepsilon_{min})\$")
ax.set_title("Methods to estimate \$\\varepsilon_{max}\$")
ax.set_yticks(0:4:20)
wsave(plotsdir("comparison", "maxε"), fig)
