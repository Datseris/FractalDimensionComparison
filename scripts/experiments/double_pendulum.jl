using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

using DelimitedFiles

name = "double_pendulum2.csv"
file = datadir("experimental", name)
M = readdlm(file, ',')

x = M[:, 5][1:15000]
y = M[:, 6][1:15000]

fig, axs = subplots(1,2; figsize = (10, 5))
axs[1].plot(x)
axs[2].plot(x,y)
axs[1].plot(y)

writedlm(datadir("experimental", "doublependulum.txt"), [x y])

estimate_embedding_plot(x)
estimate_embedding_plot(y)
# Mutual information does NOT have a clear minimum
# the value of 50 is given approximately
# Cao's method provides embedding dimension of 5

_, τ_vals_x = pecuzal_embedding(x .+ 1e-9randn(length(x)); τs = 0:100, w = 5)
_, τ_vals_y = pecuzal_embedding(y .+ 1e-9randn(length(x)); τs = 0:100, w = 5)
# Astonishingly Pecuzal method gives very similar τ for both input timeseries
τx = [0, 51, 25, 39, 12]
τx = [0, 60, 30, 46, 14]

# Let's see how it does for the dual input
A = Dataset(x .+ 1e-9randn(length(x)), y .+ 1e-9randn(length(x)))
_, τ_vals_d, ts_vals = pecuzal_embedding(A; τs = 0:100, w = 5)

# It gives a 6dimension embedding
τd = [0, 61, 31, 47, 31, 0]
jd = [2, 2, 2, 2, 1, 1]

# # PECUZAL method gives 4d embedding with τ = [0, 6, 3, 14]
#
# # %%
# writedlm(datadir("experimental", "roessler_all.txt"), M)
# writedlm(datadir("experimental", "roessler_embed.txt"), x2)
#
# # %% Okay, this will take tremendous amount of computation time,
# # using all 28 timeseries...
# # Y, τ_vals, ts_vals, = pecuzal_embedding(X1)
