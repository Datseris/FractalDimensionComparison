#=
The goal here is to provide a demonstration of the method we use to
extract the "fractal dimension", the slope of some quantity vs ε
=#

using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

N = 20_000
D = 5
dt = 0.1
ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = 24.0)
εs = ℯ .^ range(-1, 4, length = 40)
X = trajectory(ds, N*dt; dt, Ttr = 100.0)

midpoints(x) = [(x[i] + x[i+1])/2 for i in 1:length(x)-1]
localslopes(x, y) = [(y[i+1] - y[i])/(x[i+1] - x[i]) for i in 1:length(x)-1]

midpoints(x, n::Int) = [(x[i] + x[i+n])/2 for i in 1:length(x)-n]
localslopes(x, y, n::Int) = [linreg(x[i:i+n], y[i:i+n])[2] for i in 1:length(x)-n]


q = 2
H = genentropy.(Ref(X), εs; q)
C = boxed_correlationsum(X, εs; w = 10, q)

# i = findfirst(z -> z > 0, diff(-H)) - 1
# Y = -H

# i = findfirst(z -> z > 0, C)
# Y = log.(C)
# %%
fig, axs = subplots(3, 2; sharex = true, figsize = (figx, 3figy/4))
for (k, Y) in enumerate((-H, log.(C)))
    i = k==1 ? (findfirst(z -> z > 0, diff(Y)) - 1) : findfirst(isinf, Y)
    i = isnothing(i) ? 1 : max(1, i)
    x, y = log.(εs)[i:end], (Y[i:end])

    tol = 0.2
    dxi = 2

    lrs, slopes = linear_regions(x, y; tol, dxi)

    axs[1, k].plot(x, y)
    for i in 1:length(lrs)-1
        axs[1, k].plot(
            x[lrs[i]:lrs[i+1]], y[lrs[i]:lrs[i+1]];
            marker = "o", ms = 5
        )
    end

    (i1, i2), s = linear_region(x, y; tol, dxi)

    for j in 1:3
        axs[j, k].axvspan(x[i1], x[i2]; alpha = 0.25, color = "gray")
    end

    em = midpoints(x)
    δ = localslopes(x,y)
    axs[2, k].plot(em, δ, marker = "s", ms = 5)

    em = midpoints(x, 5)
    δ = localslopes(x,y, 5)
    axs[3, k].plot(em, δ, marker = "s", ms = 5)
end

axs[1,1].set_ylabel("\$-H_{$q}\$")
axs[1,1].set_yticks(0:-3:-9)
axs[1,2].set_ylabel("\$\\log(C_{$q})\$")
axs[1,2].set_yticks(0:-5:-15)
for i in 2:3, j in 1:2
    axs[i,j].set_yticks(0:2:6)
    axs[i,j].set_ylabel("\$\\delta_$(i==2 ? 1 : 5)\$")
    axs[3,j].set_xlabel("\$\\log(\\varepsilon)\$")
end

fig.tight_layout(;pad = 0.3)
# wsave(plotsdir("paper", "demonstration"), fig)
