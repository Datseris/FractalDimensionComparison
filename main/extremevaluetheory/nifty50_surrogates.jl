# Surrogate analysis of nifty50
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using Random: Xoshiro
using TimeseriesSurrogates
using DelimitedFiles
using DelayEmbeddings: genembed
import ProgressMeter

# Load experimental data
name = "nifty50"
file = datadir("experimental", name*".txt")
x = vec(readdlm(file))

# Delay times for embeddings, estimated independently
τs = [j*43 for j in 0:5]

# Parameter for EVT dimension
p = 0.98

# Surrogate method and generator
surromethod = RandomFourier()
sgen = surrogenerator(x, surromethod)

# Dimension for real data
X = genembed(x, τs)
Δx = extremevaltheory_dim(X, p)

# Surrogate dimensions
n = 100 # number of surrogates
progress = ProgressMeter.Progress(n)
Δs = map(1:n) do i
    s = sgen()
    S = genembed(s, τs)
    Δ = extremevaltheory_dim(S, p; show_progress = false)
    ProgressMeter.next!(progress)
    return Δ
end

# %% Plot
fig, ax = hist(Δs; label = "surrogates", color = (to_color(COLORS[2]), 0.75))
vlines!(ax, Δx; linestyle = :dash, color = Cycled(4), label = "original")
axislegend(ax)
ax.title = "surrogate analysis: $(name)"
ylims!(ax, (0, nothing))
display(fig)
wsave(plotsdir("paper", "evt_surrogates"), fig)
