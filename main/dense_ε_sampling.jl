using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

data = :henon_chaotic
qC = 2

data_producing_function = getfield(Data, data)
X = data_producing_function(; N = 100_000)

eC = estimate_boxsizes(X; k = 256, base = MathConstants.e)

@time C = correlationsum(X, eC; q = qC, w = 0, show_progress=true)

fig, ax = lines(log.(eC), log.(C))
ax.title = "Densely sampled ε; $(data)"
ax.xlabel = "log(ε)"
ax.ylabel = "log(C_2)"
fig