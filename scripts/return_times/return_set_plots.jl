to = Systems.towel()
tr = trajectory(to, 5000; Ttr = 10)
u0 = tr[3000]

# With these boxes, in the first 5 steps, the trajectory enters the y and z range
# but not the x range. Therefore it should NOT enter the box. See figure!
εs = [
    SVector(0.05, 0.05, 0.125),
    SVector(0.005, 0.005, 0.025),
]

# Visual guidance
using PyPlot
tr0 = trajectory(to, 5, u0)
fig, axs = subplots(1,3)
comb = ((1, 2), (1, 3), (2, 3))
for i in 1:3
    j, k = comb[i]
    ax = axs[i]
    ax.scatter(tr[:, j], tr[:, k], s = 2, color = "C$(i-1)")
    ax.scatter([u0[j]], [u0[k]], s = 20, color = "k")
    ax.plot(tr0[:, j], tr0[:, k], color = "k")
    for l in 1:length(εs)
        rect = matplotlib.patches.Rectangle(
        u0[[j, k]] .- εs[l][[j, k]], 2εs[l][j], 2εs[l][k],
        alpha = 0.25, color = "k"
        )
        ax.add_artist(rect)
    end
end


ro = Systems.roessler()
u0 = SVector(
    4.705494942754781,
    -10.221120945130545,
    0.06186563933318555,
)

εs = [
    SVector(0.05, 0.05, 0.125),
    SVector(0.005, 0.005, 0.025),
]

εs = sort!(ℯ .^ (-3:0.5:0); rev = true)

using PyPlot
tr0 = trajectory(ro, 50, u0)
fig, axs = subplots(1,3)
comb = ((1, 2), (1, 3), (2, 3))
for i in 1:3
    j, k = comb[i]
    ax = axs[i]
    ax.plot(tr[:, j], tr[:, k], lw = 0.5, color = "C$(i-1)", alpha = 0.5)
    ax.scatter([u0[j]], [u0[k]], s = 20, color = "k")
    ax.plot(tr0[:, j], tr0[:, k], color = "k", lw = 1.0, ls = "--")
    if eltype(εs[1]) <: Vector
        for l in 1:length(εs)
            rect = matplotlib.patches.Rectangle(
            u0[[j, k]] .- εs[l][[j, k]], 2εs[l][j], 2εs[l][k],
            alpha = 0.25, color = "k"
            )
            ax.add_artist(rect)
        end
    else
        for l in 1:length(εs)
            circ = matplotlib.patches.Circle(
                u0[[j, k]], εs[l]; alpha = 0.25, color = "k"
            )
            ax.add_artist(circ)
        end
    end
end
