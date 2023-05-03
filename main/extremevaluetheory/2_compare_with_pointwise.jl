using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

N = Int(1e5)
p = 0.99



# case 1: chaotic henon
data = :henon_chaotic

X = getfield(Data, data)(; N)

es = estimate_boxsizes(X)
Ds = pointwise_dimensions(X, es)

x, y = columns(X)
scatter(x, y; color = Ds, colorrange = (0.5, 1.5), markersize = 4.0)

DEs = extremevaltheory_dims_persistences(X, 0.99; compute_persistence = false)

# %%
fig, axs = axesgrid(1, 2)
