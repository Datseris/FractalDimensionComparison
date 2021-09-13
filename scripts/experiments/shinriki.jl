using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))

using DelimitedFiles
n = 2
name = "Shinriki_data_$n"
file = datadir("experimental", name*".txt")
X = readdlm(file, ','; skipstart = 21) # for Shinriki data
L = size(X)[1]
x = X[1:2:L÷10, 2]

newname = "shinriki$(n)"
writedlm(datadir("experimental", newname*".txt"), X[1:2:L÷10, 2:end])

# Run this after `x` is obtained
fig = figure()
ax = subplot(2,1,1)
ax.plot(x)

mi = selfmutualinfo(x, 0:100)
w = estimate_delay(x, "mi_min")
ax = subplot(2,2,3)
ax.plot(0:100, mi; label = "mut. inf.")
ax.axvline(w; color = "C2")
ax.legend()

Es = delay_afnn(x, w, 2:10)
ax = subplot(2,2,4)
ax.plot(2:10, Es; label = "Cao")
ax.legend()

# Y, τs = pecuzal_embedding(x .+ 1e9randn(length(x)); w, τs = 0:1:4w)
