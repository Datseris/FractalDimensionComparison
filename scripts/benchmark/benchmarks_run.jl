# Produce benchmarks comparing all functions with each other.
# Save output into a data file in data/benchmarks.jld2
using DrWatson
@quickactivate :FractalDimensionComparison
include(srcdir("style.jl"))
using BenchmarkTools, Statistics

# So that we can write a unique loop over all estimators
genentropyf(X, εs::AbstractVector) = genentropy.(Ref(X), εs)
buenooroviof(X, εs) = estimate_r0_buenoorovio(X)
takensf(X, εs) = takens_best_estimate(X, maximum(εs))
fixedmassf(X, εs) = correlationsum_fixedmass(X, length(εs))
moltenof(X, εs) = molteno_boxing(X)

Ds = 4:10
Ns = 10.0 .^ (3.3333:0.3333:5.3333) .|> round .|> Int

methods = [
    "buenoorovio r0",
    "entropy",
    "boxed large ε",
    "boxed small ε",
    "boxed small ε + r0",
    "Takens best",
    "fixed mass",
    "molteno",
]

k = 16 # amount of values in ε range
base = MathConstants.e
N_index_for_D = 6
D_index_for_N = 1

function get_times(f, X, εs)
    bm = @benchmark $f($X, $εs)
    return median(bm).time/1e9
end

function create_times(X)
    B = zeros(length(methods))
    mi, ma = minmaxima(X)
    εmax = mean(ma - mi)
    r0, ε0 = estimate_r0_buenoorovio(X)

    B[1] = get_times(buenooroviof, X, 0)

    upper_large = log(base, εmax)
    upper_small = log(base, r0)
    lower = log(base, ε0)

    small_range = float(base) .^ range(lower+1, upper_small; length = k)
    if !issorted(small_range)
        println("small range not sorted")
        small_range = float(base) .^ range(lower+1, lower+4; length = k)
    end
    large_range = float(base) .^ range(lower+1, upper_large-2; length = k)
    if !issorted(large_range)
        println("large range not sorted")
        large_range = float(base) .^ range(lower+1, upper_large+4; length = k)
    end

    B[2] = get_times(genentropyf, X, large_range)
    B[3] = get_times(boxed_correlationsum, X, large_range)
    B[4] = get_times(boxed_correlationsum, X, small_range)
    B[5] = B[4] + B[1]
    B[6] = get_times(takensf, X, large_range)
    B[7] = get_times(fixedmassf, X, large_range)
    B[8] = get_times(moltenof, X, large_range)
    return B
end

# benchmark computation versus increasing N
data = Data.lorenz96_chaotic(; D = Ds[D_index_for_N], N = maximum(Ns))
BN = zeros(length(Ns), length(methods))
for (i, N) in enumerate(Ns)
    @show i, N
    X = data[1:N]
    BN[i, :] = create_times(X)
end

# benchmark computation versus increasing D
BD = zeros(length(Ds), length(methods))
for (i, D) in enumerate(Ds)
    @show i, D
 X = Data.lorenz96_chaotic(; D, N = Ns[N_index_for_D])
    BD[i, :] = create_times(X)
end

# sname = savename("benchmarks", @dict(M, S, T, k), "jld2")

wsave(datadir("benchmarks", "benchmarks.jld2"), (@strdict Ds Ns methods BN BD D_index_for_N N_index_for_D k base))
