#=
This script generates the "pitfall" data from raw input data
=#

using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))
using DelimitedFiles

# %% nifty 50
name = "nifty50.csv"
file = datadir("pitfall", name)
X = readdlm(file, ','; skipstart = 1)
xx = X[:, 2]
# fix "null" nonsense
for (i, e) in enumerate(xx)
    e isa AbstractFloat && continue
    js = i+1:length(xx)
    j = findfirst(ι -> xx[ι] isa AbstractFloat, js)
    xx[i] = xx[js[j]]
end
x = float.(xx)
file = datadir("pitfall", name)
writedlm(datadir("experimental", "nifty50.txt"), x)

# Run this after `x` is obtained
# Should we de-trend here?
estimate_embedding(x, "nifty50")


# %% Vostok
name = "vostok.csv"
file = datadir("pitfall", name)
X = readdlm(file, '\t')
tim = X[:, 1]
δT = X[:, 2]

using Dierckx
# interpolation object:
spl = Spline1D(tim, δT; k=3, bc="extrapolate")
# equi-spaced time vector:
t = range(minimum(tim), maximum(tim); length = 2000)
T = spl(t)
writedlm(datadir("experimental", "vostok.txt"), T)

# Run this after `x` is obtained
# Should we de-trend here?
estimate_embedding(T, "vostok")
