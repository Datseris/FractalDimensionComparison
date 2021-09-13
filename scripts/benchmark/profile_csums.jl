using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))
using BenchmarkTools

X = Data.henon_chaotic(N = 10_000)
es = estimate_boxsizes(X)
qC = 2

for qC in (2, 3, 4)
    @show qC
    bmark = @benchmark correlationsum($X, $es; q = qC)
    bmark_boxed = @benchmark boxed_correlationsum($X, $es; q = qC)
    
    @show median(bmark).time
    @show median(bmark_boxed).time
    @show median(bmark).allocs
    @show median(bmark_boxed).allocs
    println()
end
