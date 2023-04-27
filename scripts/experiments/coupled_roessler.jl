using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

using DelimitedFiles

name = "R1_ST_10_1.dat"

file = datadir("exp_pro", name)
M = readdlm(file)


X1 = Dataset(M) # Every single oscillator as one dimension

x2 = M[:, 1]
estimate_embedding_plot(x2)
# From the above we obtain embedding dimension of ~6 and delay time of 7

Y, τ_vals = pecuzal_embedding(x2; τs = 0:50, w = 6)
# PECUZAL method gives 4d embedding with τ = [0, 6, 3, 14]

# %%
writedlm(datadir("experimental", "roessler_all.txt"), M)
writedlm(datadir("experimental", "roessler_embed.txt"), x2)

# %% Okay, this will take tremendous amount of computation time,
# using all 28 timeseries...
# Y, τ_vals, ts_vals, = pecuzal_embedding(X1)
