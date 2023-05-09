using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

# Plot state space attractor colorcoded
data = :lorenz63_chaotic

X = getfield(Data, data)(; N, Δt)

es = estimate_boxsizes(X)
Ds = pointwise_dimensions(X, es)


DEs, _ = extremevaltheory_dims_persistences(X, 0.99; compute_persistence = false)

# %%
x, y = X[:, 1], X[:, 3]

fig, axs = axesgrid(1, 2)

markersize = 5
cmap = :viridis
cmapkw = (colormap = cmap, lowclip = :orange, highclip = :red, colorrange = (1.0, 3.0))

scatter!(axs[1], x, y; markersize, color = Ds, cmapkw...)
scatter!(axs[2], x, y; markersize, color = DEs, cmapkw...)

cbar = fig[0, :] = Colorbar(fig; cmapkw..., vertical = false)

# Distribution difference as violin plot
axd = Axis(fig[2, :])
rowsize!(fig.layout, 2, Relative(1/3))

density!(axd, Ds; color = (to_color(COLORS[1]), 0.5), strokewidth = 3, strokecolor = COLORS[1])
density!(axd, DEs; color = (to_color(COLORS[2]), 0.5), strokewidth = 3, strokecolor = COLORS[2])

fig