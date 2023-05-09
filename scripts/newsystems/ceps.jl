using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using DelimitedFiles, DelayEmbeddings

datapath = datadir("ceps")

files = sort!(readdir(datapath))

file = files[1]

X = vec(readdlm(joinpath(datapath, file), ','))
X = X[1:5:end]

fig = scatterlines(X; axis = (title = file[1:14],), markersize = 10)
display(fig)

𝒟, τ, E = optimal_traditional_de(X)
@show dimension(𝒟)
@show τ

Δc = grassberger_proccacia_dim(𝒟; w = τ)
@show Δc
Δh = generalized_dim(𝒟)
@show Δh

# Y, tau_vals, = pecuzal_embedding(X; max_cycles = 6, w = 10, verbose = false)
# @show tau_vals

# Y, τ_vals, ts_vals, = pecuzal_embedding(x)
