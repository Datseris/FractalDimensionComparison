#=
The goal here is to provide a demonstration of the method we use to
extract the "fractal dimension", the slope of some quantity vs ε
=#

using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using PredefinedDynamicalSystems
N = 20_000
D = 5

X = Data.lorenz96_chaotic(; D = 5, N = 20_000)
εs = ℯ .^ range(-4, 3, length = 40)

midpoints(x) = [(x[i] + x[i+1])/2 for i in 1:length(x)-1]
localslopes(x, y) = [(y[i+1] - y[i])/(x[i+1] - x[i]) for i in 1:length(x)-1]

midpoints(x, n::Int) = [(x[i] + x[i+n])/2 for i in 1:length(x)-n]
localslopes(x, y, n::Int) = [linreg(x[i:i+n], y[i:i+n])[2] for i in 1:length(x)-n]

q = 2
H = [entropy(Renyi(q, MathConstants.e), ValueHistogram(ε), X) for ε in εs]
C = boxed_correlationsum(X, εs; w = 10, q)

# i = findfirst(z -> z > 0, diff(-H)) - 1
# Y = -H

# i = findfirst(z -> z > 0, C)
# Y = log.(C)
# %%
fig, axs = axesgrid(3, 2; sharex = true)

for (k, Y) in enumerate((-H, log.(C)))
    i = k==1 ? (findfirst(z -> z > 0, diff(Y)) - 1) : findfirst(isinf, Y)
    i = isnothing(i) ? 1 : max(1, i)
    x, y = log.(εs)[i:end], (Y[i:end])

    tol = 0.2
    dxi = 2

    lrs, slopes = linear_regions(x, y; tol, dxi)

    lines!(axs[1, k], x, y)
    for i in eachindex(lrs)
        scatterlines!(axs[1, k],
            x[lrs[i]], y[lrs[i]];
            markersize = 10, color = Cycled(i),
        )
    end

    reg, s = linear_region(x, y; tol, dxi)

    for j in 1:3
        vspan!(axs[j, k], x[reg[1]], x[reg[end]]; color = ("gray", 0.25))
    end

    em = midpoints(x)
    δ = localslopes(x,y)
    scatterlines!(axs[2, k], em, δ; markersize = 10)

    em = midpoints(x, 5)
    δ = localslopes(x,y, 5)
    scatterlines!(axs[3, k], em, δ; markersize = 10)
end

axs[1,1].ylabel = L"-H_{%$q}"
axs[1,1].yticks = 0:-3:-9
axs[1,2].ylabel = L"\log(C_{%$q})"
axs[1,2].yticks = 0:-5:-15
for j in 1:2
    axs[3,j].xlabel = L"\log(\varepsilon)"
    for i in 2:3
        axs[i,j].yticks = 0:2:6
        axs[i,j].ylabel = L"\delta_%$(i==2 ? 1 : 5)"
    end
end

fig
wsave(plotsdir("paper", "demonstration"), fig)
