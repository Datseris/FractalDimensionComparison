# Produce benchmarks comparing all functions with each other
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))
using BenchmarkTools, Statistics

# So that we can write a unique loop over all estimators
genentropyf(X, εs::AbstractVector) = genentropy.(Ref(X), εs)
buenooroviof(X, εs) = estimate_r0_buenoorovio(X)
takensf(X, εs) = takens_best_estimate(X, maximum(εs))

Ds = 4:10 #10 # Dimensions
Ns = 10.0 .^ (3.3333:0.3333:5.3333) .|> round .|> Int

methods = [
    "estimate r0",
    "entropy",
    "boxed large ε",
    "boxed small ε",
    "boxed small ε + r0",
    "Takens best",
]

k = 20 # amount of values in ε range
base = MathConstants.e
Ni = 6
Di = 1

function get_times_b(f, X, εs)
    bm = @benchmark f($X, $εs)
    return median(bm).time/1e9
end

function get_times_e(f, X, εs)
    bm = @elapsed f(X, εs)
    return bm
end

get_times = get_times_b

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
        println("small range not sorted")
        large_range = float(base) .^ range(lower+1, upper_large+4; length = k)
    end

    B[2] = get_times(genentropyf, X, large_range)
    B[3] = get_times(boxed_correlationsum, X, large_range)
    B[4] = get_times(boxed_correlationsum, X, small_range)
    B[5] = B[4] + B[1]
    B[6] = get_times(takensf, X, large_range)
    return B
end

data = Data.lorenz96_chaotic(; D = Ds[Di], N = maximum(Ns))
BN = zeros(length(Ns), length(methods))
for (i, N) in enumerate(Ns)
    @show i, N
    X = data[1:N]
    BN[i, :] = create_times(X)
end

BD = zeros(length(Ds), length(methods))
for (i, D) in enumerate(Ds)
    @show i, D
    X = Data.lorenz96_chaotic(; D, N = Ns[Ni])
    BD[i, :] = create_times(X)
end

# sname = savename("benchmarks", @dict(M, S, T, k), "jld2")

wsave(datadir("benchmarks", "benchmarks.jld2"), (@strdict Ds Ns methods BN BD Di Ni k base))
