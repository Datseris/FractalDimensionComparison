using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, BenchmarkTools

import ChaosTools: float_to_int, molteno_boxing, _molteno_boxing, estimate_r0_buenoorovio, correlation_boxing

function float_to_int(data::Dataset{D,T}) where {D, T}
    N = length(data)
    mins, maxs = minmaxima(data)
    sizes = maxs .- mins
    ε0 = maximum(sizes)
    I = T(typemax(UInt64))
    # Let f:[min,max] -> [0+eps(T),1-eps(T)]*typemax(UInt64), then f(x) = m*x + b
    m = (1 - 2eps(T)) * (I - 2eps(I)) / ε0
    b = eps(I) .- mins .* m

    res = Vector{SVector{D,UInt64}}()
    sizehint!(res, N)
    for x in data
        int_val = floor.(UInt64, m .* x + b)
        push!(res, int_val)
    end
    Dataset(res), ε0
end

return_T_eps(X::Dataset{D,T}) where {D,T} = eps(T(typemax(UInt64)))

function molteno_boxing(data::Dataset, k0 = 10)
    integers, ε0 = float_to_int(data)
    boxes = _molteno_boxing(integers, k0)
    εs = ε0 ./ 2 .^ (1:length(boxes))
    return boxes, εs
end

function estimate_r0_buenoorovio(data)
    r0 = zero(eltype(data))
    mini, maxi = minmaxima(data)
    N = length(data)
    m = DelayEmbeddings.dimension(data)
    R = maximum(maxi .- mini)
    # The possibility of a bad pick exists, if so, the calculation is repeated.
    ν = zero(eltype(data))
    # Sample N/10 datapoints out of data for rough estimate of effective size.
    data_sample2 = data[unique(rand(1:N, N÷10))] |> Dataset
    r_ℓ = R / 10
    η_ℓ = length(correlation_boxing(data_sample2, r_ℓ)[1])
    r0 = zero(eltype(data))
    while true
        # Sample √N datapoints for rough dimension estimate
        data_sample1 = data[unique(rand(1:N, ceil(Int, sqrt(N))))] |> Dataset
        # Define logarithmic series of radii.
        lower = log10(min_pairwise_distance(data_sample1)[2])
        εs = 10 .^ range(lower, stop = log10(R), length = 16)
        # Estimate ν from a sample using the Grassberger Procaccia algorithm.
        cm = correlationsum(data_sample1, εs)
        ν = linear_region(log.(εs), log.(cm), tol = 0.5)[2]
        # Estimate the effictive size of the chaotic attractor.
        ℓ = r_ℓ * η_ℓ^(1/ν)
        # Calculate the optimal number of filled boxes according to Bueno-Orovio
        η_opt = N^(2/3) * ((3^ν - 1/2) / (3^m - 1))^(1/2)
        # The optimal box size is the effictive size divided by the box number # to the power of the inverse dimension.
        r0 = ℓ / η_opt^(1/ν)
        !isnan(r0) && break
    end
    r0
end

@btime()

# %% TODO: Also compare an embedded set for different `d` and `τ`.

# %% sensitivity to data dimensionality
N, dt = 2000, 0.1
fig, axs = subplots(4,1; sharex = true, figsize = (10,16))
for (j, D) in enumerate(4:3:16)
    color = "C$(j-1)"
    ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = 8.0)
    X = trajectory(ds, N*dt; dt = dt, Ttr = 100.0)
    r0 = estimate_r0_buenoorovio(X)

    εs = 10.0 .^ range(-4 + D/4, 2, length = 32)
    bεs = r0 .* 10 .^ range(log10(min_pairwise_distance(X)[2]/r0), stop = 0, length = 16)
    m, mεs = molteno_boxing(X)

    H2 = genentropy.(2, εs, Ref(X); base = 10.0)
    H1 = genentropy.(2, m, base = 10.0)
    C1 = correlationsum(X, εs; w = 5)
    C2 = boxed_correlationsum(X, bεs, r0)
    axs[1].plot(log10.(εs), -H2, label = "\$ D = $D\$", color = color)
    if D == 4 || D == 7
        i = findfirst(z -> z > 0, 10 .^ (-H1))
        x, y = log10.(mεs)[i:end], -H1[i:end]
        is, d = linear_region(x, y)
        axs[2].plot(x, y, zorder = 1, color = color, label = "\$ d = $(round(d, digits = 2)) \$")
        axs[2].plot([x[is[1]], x[is[end]]], [y[is[1]], y[is[end]]],
        zorder = 2, ms = 5, marker = "o", ls = "None", color = color)
    else
        axs[2].plot(log10.(mεs), -H1, color = color)
    end
    axs[4].plot(log10.(bεs), log10.(C2))
    axs[4].set_xlabel("\$\\log_{10}(\\epsilon)\$")
    axs[4].set_ylabel("\$\\log_{10}(C_2)\$")
    axs[3].set_ylabel("\$\\log_{10}(C_1)\$")
    axs[2].set_ylabel("\$-H_1\$")
    axs[1].set_ylabel("\$-H_2\$")
    for (number, (set, es)) in enumerate(zip([C1, C2], [εs, bεs]))
        i = findfirst(z -> z > 0, set)
        x, y = log10.(es)[i:end], log10.(set)[i:end]
        is, d = linear_region(x, y)
        axs[2+number].plot(x, y, zorder = 1, color = color, label = "\$ d = $(round(d, digits = 2)) \$")
        axs[2+number].plot([x[is[1]], x[is[end]]], [y[is[1]], y[is[end]]],
        zorder = 2, ms = 5, marker = "o", ls = "None", color = color)
    end
end
for ax in axs; ax.legend(loc = "upper left"); end
fig.tight_layout()
axs[1].tick_params(labelbottom=false)
axs[1].set_title("increasing data dimension (Lorenz-96, N=$(N))")
wsave("plots/comparing_datadim.png", fig)
# comment: N = 100_000 took so much time on my pc I had to shut it down.


# %% Sensititivy to initial condition

N, dt, D = 5000, 0.1, 4
fig, axs = subplots(4,1; sharex = true, figsize = (10,16))
for (j, __D) in enumerate(4:3:16)
    color = "C$(j-1)"
    ds = Systems.lorenz96(D; F = 8.0)
    X = trajectory(ds, N*dt; dt = dt, Ttr = 10.0)

    r0 = estimate_r0_buenoorovio(X)

    εs = 10.0 .^ range(-4 + D/4, 2, length = 32)
    bεs = r0 .* 10 .^ range(log10(min_pairwise_distance(X)[2]/r0), stop = 0, length = 16)
    m, mεs = molteno_boxing(X)

    H2 = genentropy.(2, εs, Ref(X); base = 10.0)
    H1 = genentropy.(2, m, base = 10.0)
    C1 = correlationsum(X, εs; w = 5)
    C2 = boxed_correlationsum(X, bεs, r0)
    axs[1].plot(log10.(εs), -H2, label = nothing, color = color)
    axs[2].plot(log10.(mεs), -H1, color = color)
    axs[4].plot(log10.(bεs), log10.(C2))
    axs[4].set_xlabel("\$\\log_{10}(\\epsilon)\$")
    axs[4].set_ylabel("\$\\log_{10}(C_2)\$")
    axs[3].set_ylabel("\$\\log_{10}(C_1)\$")
    axs[2].set_ylabel("\$-H_1\$")
    axs[1].set_ylabel("\$-H_2\$")
    for (number, (set, es)) in enumerate(zip([10 .^ -H1, C1, C2], [mεs, εs, bεs]))
        i = findfirst(z -> z > 0, set)
        x, y = log10.(es)[i:end], log10.(set)[i:end]
        is, d = linear_region(x, y)
        axs[1+number].plot(x, y, zorder = 1, color = color, label = "\$ d = $(round(d, digits = 2))\$")
        axs[1+number].plot([x[is[1]], x[is[end]]], [y[is[1]], y[is[end]]],
        zorder = 2, ms = 5, marker = "o", ls = "None", color = color)
    end
end
for ax in axs; ax.legend(loc = "upper left"); ax.grid(); end
axs[1].tick_params(labelbottom=false)
axs[1].set_title("different initial conditions (Lorenz-96, N=$(N))")
fig.savefig("plots/comparing_initcon.png")


# %% Sensititivy to trajectory length

dt = 0.1
fig, axs = subplots(4,1; sharex = true, figsize = (10, 16))
for (j, N) in enumerate([500, 1000, 5000, 10000, 20000, 50000])
    color = "C$(j-1)"
    ds = Systems.lorenz()
    X = trajectory(ds, N*dt; dt = dt, Ttr = 10.0)

    r0 = estimate_r0_buenoorovio(X)

    εs = 10.0 .^ range(-4 + 3/4, 2, length = 40)
    bεs = r0 .* 10 .^ range(log10(min_pairwise_distance(X)[2]/r0), stop = 0, length = 40)
    m, mεs = molteno_boxing(X)

    H2 = genentropy.(2, εs, Ref(X); base = 10.0)
    H1 = genentropy.(2, m, base = 10.0)
    C1 = [correlationsum(X, ε; w = 5) for ε in εs]
    C2 = boxed_correlationsum(X, bεs, r0)
    axs[1].plot(log10.(εs), -H2, label = "\$N = $N\$", color = color)
    #axs[2].plot(log10.(mεs), -H1, color = color)
    #axs[4].plot(log10.(bεs), log10.(C2))
    axs[4].set_xlabel("\$\\log_{10}(\\epsilon)\$")
    axs[4].set_ylabel("\$\\log_{10}(C_2)\$")
    axs[3].set_ylabel("\$\\log_{10}(C_1)\$")
    axs[2].set_ylabel("\$-H_1\$")
    axs[1].set_ylabel("\$-H_2\$")
    for (number, (set, es)) in enumerate(zip([10 .^ -H1, C1, C2], [mεs, εs, bεs]))
        i = findfirst(z -> z > 0, set)
        x, y = log10.(es)[i:end], log10.(set)[i:end]
        is, d = linear_region(x, y)
        axs[1+number].plot(x, y, zorder = 1, color = color, label = "\$ d = $(round(d,digits=2))\$")
        axs[1+number].plot([x[is[1]], x[is[end]]], [y[is[1]], y[is[end]]],
        zorder = 2, ms = 5, marker = "o", ls = "None", color = color)
    end
end
for ax in axs; ax.legend(loc = "upper left"); ax.grid(); end
axs[1].tick_params(labelbottom=false)
axs[1].set_title("different trajectory lengths (Lorenz63)")
fig.savefig("plots/comparing_length.png")


# %% sensitivity to underlying chaotic set
N, dt = 10000, 0.1
D = 7

lo67 = Systems.lorenz()
X1 = trajectory(lo67, N*dt; dt = dt, Ttr = 10.0)

lo96 = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = 8.0)
X2 = trajectory(lo96, N*dt; dt = dt, Ttr = 10.0)

he = Systems.henon()
X3 = trajectory(he, N; Ttr = 10)

to = Systems.towel()
X4 = trajectory(to, N; Ttr = 10)

sm = Systems.coupledstandardmaps(4; ks = fill(0.5, 4), Γ = 0.5)
X5 = trajectory(sm, N; Ttr = 100)

sm = Systems.coupledstandardmaps(4; ks = fill(0.5, 4), Γ = 0.5)
X5 = trajectory(sm, N; Ttr = 100)

using Distributions # for multivariate gaussian
dist = MvNormal(6, 0.1)
X6 = Dataset(rand(dist, N)')

names = ["Lorenz63", "Lorenz96 (7)", "Hénon", "towel", "coupled SM (8)", "gaussian (6)"]
lss = ["-", "--", ":", "-."]
fig, axs = subplots(4,1; sharex = true, figsize = (10,16))

for (j, X) in enumerate((X1, X2, X3, X4, X5, X6))
    color = "C$(j-1)"
    ls = lss[mod1(j, 4)]
    name = names[j]

    r0 = estimate_r0_buenoorovio(X)

    εs = 10.0 .^ range(-4 + 3/4, 2, length = 40)
    bεs = r0 .* 10 .^ range(log10(min_pairwise_distance(X)[2]/r0), stop = 0, length = 40)
    m, mεs = molteno_boxing(X)

    H2 = genentropy.(2, εs, Ref(X); base = 10.0)
    H1 = genentropy.(2, m, base = 10.0)
    C1 = [correlationsum(X, ε; w = 5) for ε in εs]
    C2 = boxed_correlationsum(X, bεs, r0)
    axs[1].plot(log10.(εs), -H2, label = name, color = color, ls = ls)
    #axs[2].plot(log10.(mεs), -H1, color = color)
    #axs[4].plot(log10.(bεs), log10.(C2))
    axs[4].set_xlabel("\$\\log_{10}(\\epsilon)\$")
    axs[4].set_ylabel("\$\\log_{10}(C_2)\$")
    axs[3].set_ylabel("\$\\log_{10}(C_1)\$")
    axs[2].set_ylabel("\$-H_1\$")
    axs[1].set_ylabel("\$-H_2\$")
    for (number, (set, es)) in enumerate(zip([10 .^ -H1, C1, C2], [mεs, εs, bεs]))
        i = findfirst(z -> z > 0, set)
        x, y = log10.(es)[i:end], log10.(set)[i:end]
        is, d = linear_region(x, y)
        axs[1+number].plot(x, y, zorder = 1, color = color, label = "\$ d = $(round(d, digits = 2))\$", ls = ls)
        axs[1+number].plot([x[is[1]], x[is[end]]], [y[is[1]], y[is[end]]],
        zorder = 2, ms = 5, marker = "o", ls = "None", color = color)
    end
end
for ax in axs; ax.legend(loc = "upper left"); ax.grid(); end
axs[1].tick_params(labelbottom=false)
axs[1].set_title("different datasets (N=$(N))")
fig.savefig("plots/comparing_differentedata.png")

# %% computation times for different lengths

dt = 0.1
Ns = 500:500:10000
times_lo67 = [[] for i in 1:4]
times_he = [[] for i in 1:4]
for (j, N) in enumerate(Ns)
    X1 = trajectory(lo67, N*dt; dt = dt, Ttr = 10.0)
    X2 = trajectory(he, N; Ttr = 10)
    # standard generalized entropy algorithm
    t_lo67 = @benchmark begin
        ε0 = min_pairwise_distance(X1)[2]
        ε1 = minimum(maxima(X1) - minima(X1))
        εs = 10.0 .^ range(log10(ε0), log10(ε1), length = 12)
        H = genentropy.(2, εs, Ref(X1); base = 10.0)
    end
    t_he = @benchmark begin
        ε0 = min_pairwise_distance(X2)[2]
        ε1 = minimum(maxima(X2) - minima(X2))
        εs = 10.0 .^ range(log10(ε0), log10(ε1), length = 12)
        H = genentropy.(2, εs, Ref(X2); base = 10.0)
    end
    push!(times_lo67[1], minimum(t_lo67.times))
    push!(times_he[1], minimum(t_he.times))

    # molteno's generalized entropy algorithm
    t_lo67 = @benchmark begin
        boxes, εs = molteno_boxing(X1)
        H = genentropy.(2, boxes; base = 10.0)
    end
    t_he = @benchmark begin
        boxes, εs = molteno_boxing(X2)
        H = genentropy.(2, boxes; base = 10.0)
    end
    push!(times_lo67[2], minimum(t_lo67.times))
    push!(times_he[2], minimum(t_he.times))

    # classic correlation algorithm
    t_lo67 = @benchmark begin
        ε0 = min_pairwise_distance(X1)[2]
        ε1 = minimum(maxima(X1) - minima(X1))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        C = log10.(correlationsum(X1, εs; w = 5))
    end
    t_he = @benchmark begin
        ε0 = min_pairwise_distance(X2)[2]
        ε1 = minimum(maxima(X2) - minima(X2))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        C = log10.(correlationsum(X2, εs; w = 5))
    end
    push!(times_lo67[3], minimum(t_lo67.times))
    push!(times_he[3], minimum(t_he.times))

    # boxed correlation algorithm
    t_lo67 = @benchmark begin
        r0 = estimate_r0_buenoorovio(X1)
        ε0 = min_pairwise_distance(X1)[2]
        εs = 10.0 .^ range(log10(ε0), log10(r0), length = 12)
        C = log10.(boxed_correlationsum(X1, εs, r0))
    end
    t_he = @benchmark begin
        r0 = estimate_r0_buenoorovio(X2)
        ε0 = min_pairwise_distance(X2)[2]
        εs = 10.0 .^ range(log10(ε0), log10(r0), length = 12)
        C = log10.(boxed_correlationsum(X2, εs, r0))
    end
    push!(times_lo67[4], minimum(t_lo67.times))
    push!(times_he[4], minimum(t_he.times))
end
fig, axs = subplots(2,1; sharex = true, figsize = (10, 16))
names = ["Generalized dimension" "Molteno's method" "Correlation dimension" "Boxed correlation dimension"]
for i in 1:4
    color = "C$(i-1)"
    name = names[i]
    axs[1].plot(Ns, times_lo67[i], label = name, color = color)
    axs[2].plot(Ns, times_he[i], color = color)
end
axs[2].set_xlabel("\$N\$")
axs[2].set_ylabel("\$t_{Hénon}\$")
axs[1].set_ylabel("\$t_{Lorenz63}\$")


for ax in axs; ax.legend(loc = "upper left"); ax.grid(); end
axs[1].tick_params(labelbottom=false)
axs[1].set_title("minimum computation time")
fig.savefig("plots/comparing_btime.png")

# %% sensitivity to noise
