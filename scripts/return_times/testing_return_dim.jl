using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using DynamicalSystems, PyPlot, Statistics
include(srcdir("style.jl"))

# Plot dependence on time and i.c.
function plot_return_dim(ds, εs, ic, e; Tmax = 20_000_000, ici=1)
    fig, axs = subplots(2, 1; sharex = true)
    mcts = []
    for (i, T) ∈ enumerate(ic)
        println("ic $i ...")
        u0 = trajectory(ds, T; dt = T, Ttr = 100)[end]
        n, c = mean_return_times(ds, u0, εs, Tmax)
        j = findfirst(isequal(0), c)
        j = isnothing(j) ? length(c) : j-1
        push!(mcts, n[1:j])
        o, d = linreg(log.(εs[1:j]), log.(n[1:j]))
        axs[1].plot(log.(εs), log.(n);
                    marker = "o", ls="-", lw = 1.0, ms = 5, color = "C$(i)",
                    label = "i.c. $(i), \$d = $(rdspl(abs(d)))\$"
        )
    end
    for (i, T) ∈ enumerate(e)
        println("Time $i ...")
        u0 = trajectory(ds, ic[ici]; Ttr = 100)[end]
        n, c = mean_return_times(ds, u0, εs, 10^T)
        j = findfirst(isequal(0), c)
        j = isnothing(j) ? length(c) : j-1
        o, d = linreg(log.(εs[1:j]), log.(n[1:j]))
        axs[2].plot(log.(εs), log.(n);
                    marker = "o", ls="-", lw = 1.0, ms = 5, color = lighten_color("C0", 1.2i),
                    label = "\$T=10^$(T), d = $(rdspl(abs(d)))\$"
        )
    end

    minlength = minimum(length.(mcts))
    average_mct = [mean(mcts[j][i] for j in 1:length(ic)) for i in 1:minlength]
    x, y = log.(εs[1:minlength]), log.(average_mct)
    _, d = linreg(x, -y)
    # _, d = linear_region(x, -y) # this doesn't work well. The line is too jagged.
    axs[1].plot(x, y; ls = "dashed", label = "avg., \$d=$(rdspl(d))\$", lw = 3.0)

    axs[1].legend(); axs[2].legend()
    axs[1].set_ylabel("log(⟨t⟩)")
    axs[2].set_ylabel("log(⟨t⟩)")
    axs[2].set_xlabel("log(ε)")
    fig.tight_layout()
    return fig, axs
end

# %% Henon map
# first plot random initial conditions
ds = Systems.henon()
ic = (500, 1000, 1500, 2000, 2500)
e = (3, 5, 7)
εs = sort!(ℯ .^ (-9:0.5:-3); rev = true)
fig, axs = plot_return_dim(ds, εs, ic, e)
axs[1].set_title("Henon map")

tr = trajectory(ds, 50000; Ttr = 100)
figure()
scatter(columns(tr)..., s = 2)
for (i, T) ∈ enumerate(ic)
    u0 = trajectory(ds, T; Ttr = 100)[end]
    scatter(u0[1], u0[2], color = "C$i")
end

# wsave(plotsdir("returntime", "henon"), fig)

# %% Towel map
# first plot random initial conditions
ds = Systems.towel()
εs = sort!(ℯ .^ (-7:0.5:-3); rev = true)
ic = (500, 1000, 1500)
e = (4, 6, 8)

tr = trajectory(ds, 50000; Ttr = 100)
using3D(); figure();
scatter3D(columns(tr)...; s = 2)
for (i, T) ∈ enumerate(ic)
    u0 = trajectory(ds, T; Ttr = 100)[end]
    scatter3D([u0[1]], [u0[2]], [u0[3]]; color = "C$i", s = 1000)
end

fig, axs = plot_return_dim(ds, εs, ic, e)
axs[1].set_title("Towel map")
# wsave(plotsdir("returntime", "towel"), fig)

# %% Roessler system
# first plot random initial conditions
ds = Systems.roessler([1,1,1.0])
εs = sort!(ℯ .^ (-4:0.5:-1); rev = true)
ic = (563.0, 1300.0, 1000.0, 1500.0, 222.0)
e = (4, 5, 6)

tr = trajectory(ds, 1000; Ttr = 100)
using3D(); figure();
plot3D(columns(tr)...; lw = 0.5)
for (i, T) ∈ enumerate(ic)
    u0 = trajectory(ds, T; Ttr = 100)[end]
    scatter3D([u0[1]], [u0[2]], [u0[3]]; color = "C$i", s = 1000)
end

# fig, axs = plot_return_dim(ds, εs, ic, e; Tmax = 10^6, ici=2)
# axs[1].set_title("Roessler system")
# wsave(plotsdir("returntime", "roessler"), fig)

# %% Roessler system
# first plot random initial conditions
ds = Systems.henonheiles([0.0, -0.25, 0.42081, 0.0])
εs = sort!(ℯ .^ (-4:0.5:-1); rev = true)
ic = (563.0, 1300.0, 1002.0, 1500.5, 222.0)
e = (4, 5, 6)

tr = trajectory(ds, 1000; Ttr = 100)
using3D(); figure();
x, y, z = columns(tr)
plot3D(x, y, z; lw = 0.5)
for (i, T) ∈ enumerate(ic)
    u0 = trajectory(ds, T; Ttr = 100)[end]
    scatter3D([u0[1]], [u0[2]], [u0[3]]; color = "C$i", s = 1000)
end

# fig, axs = plot_return_dim(ds, εs, ic, e; Tmax = 10^6, ici=2)
# axs[1].set_title("Roessler system")
# wsave(plotsdir("returntime", "roessler"), fig)

# %% Lorenz system
# first plot random initial conditions
ds = Systems.lorenz([1,10.0,1.0])
εs = sort!(ℯ .^ (-5:0.5:-0); rev = true)
ic = (563, 1000, 1500.5)
e = (4, 5, 6)

tr = trajectory(ds, 500; Ttr = 100, dt = 0.02)
using3D(); figure();
plot3D(columns(tr)...; lw = 1.0)
for (i, T) ∈ enumerate(ic)
    u0 = trajectory(ds, T; Ttr = 100)[end]
    scatter3D([u0[1]], [u0[2]], [u0[3]]; color = "C$i", s = 1000)
end

fig, axs = plot_return_dim(ds, εs, ic, e; Tmax = 10^5, ici=2)
axs[1].set_title("Lorenz system")
# wsave(plotsdir("returntime", "lorenz"), fig)

# %% Lorenz-96 system
D = 4
ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = 24.0)
εs = sort!(ℯ .^ (-5:0.5:-0); rev = true)
ic = (26.0, 563.0, 876.8, 1000.0)
e = (4, 5, 6)

# TODO: Must increase max return time.

tr = trajectory(ds, 50; Ttr = 100.0, dt = 0.01)
using3D(); figure();
x, y, z = columns(tr)
plot3D(x, y, z; lw = 0.5)
for (i, T) ∈ enumerate(ic)
    u0 = trajectory(ds, T; Ttr = 100)[end]
    scatter3D([u0[1]], [u0[2]], [u0[3]]; color = "C$i", s = 1000)
end

fig, axs = plot_return_dim(ds, εs, ic, e; Tmax = 10^5, ici=2)
axs[1].set_title("Lorenz96, D = $D")
wsave(plotsdir("returntime", "lorenz96_D=$D"), fig)
