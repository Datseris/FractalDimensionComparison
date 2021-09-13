# Produce benchmarks comparing all functions with each other
using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))

file = wload(datadir("benchmarks", "benchmarks.jld2"))
@unpack Di, Ni, Ds, Ns, BN, BD, methods = file
BD = log10.(BD) # better to display timings in log scale
BN = log10.(BN) # better to display timings in log scale

fig, axs = subplots(2,1)

for (i, m) in enumerate(methods)
    axs[1].plot(log10.(Ns), BN[:, i]; 
    label = m, ls = LINESTYLES[i], alpha = 0.9, marker = MARKERS[i], ms = 8)
    axs[2].plot(Ds, BD[:, i], ls = LINESTYLES[i], alpha = 0.9, marker = MARKERS[i], ms = 8)
end
leg = axs[1].legend(
    bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
    ncol=3, mode="expand", borderaxespad=0,
)

axs[1].set_xlabel("\$\\log_{10}(N)\$")
axs[1].set_ylabel("\$\\log_{10}(t \\mathrm{[sec]})\$")
axs[2].set_xlabel("\$D\$")
axs[2].set_ylabel("\$\\log_{10}(t \\mathrm{[sec]})\$")
fig.tight_layout(pad=0.3)