using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using DynamicalSystems, PyPlot, Statistics
include(srcdir("style.jl"))

using OrdinaryDiffEq, DiffEqBase


# This function over-writes the currently unexported version of ChaosTools, which
# has nested `εs` but doesn't work because of diffeq error, see here:

function ChaosTools.transit_time_statistics(ds::ContinuousDynamicalSystem, u0, εs, T;
        alg = Tsit5(), diffeq...
    )

    eT = eltype(ds.t0)
    ChaosTools.check_εs_sorting(εs, length(u0))
    exits = [eT[] for _ in 1:length(εs)]
    entries = [eT[] for _ in 1:length(εs)]

    # Make the magic callback:
    crossing(u, t, integ) = ChaosTools.εdistance(u, u0, εs[1])
    negative_affect!(integ) = push!(entries[1], integ.t)
    positive_affect!(integ) = push!(exits[1], integ.t)
    cb = ContinuousCallback(crossing, positive_affect!, negative_affect!;
        save_positions = (false, false)
    )

    prob = ODEProblem(ds, (eT(0), eT(T)); u0 = u0)
    sol = solve(prob, alg;
        callback=cb, save_everystep = false, dense = false,
        save_start=false, save_end = false, diffeq...
    )
    return exits, entries
end


ro = Systems.roessler()
u0 = SVector(
    4.705494942754781,
    -10.221120945130545,
    0.06186563933318555,
)

exponents = (-4:0.5:-1)
n = Float64[]
εs = Float64[]

for e in exponents
    x = [ℯ^e]
    exits, entries = transit_time_statistics(ro, u0, x, 100.0^(-e))
    transits, returns = transit_return(exits, entries)
    m = mean.(returns)
    push!(εs, x[1])
    push!(n, m[1])
end


# plot
fig, ax = subplots()

o, d = linreg(log.(εs), log.(n))
ax.plot(log.(εs), log.(n), marker = "o", ls="-", lw = 1.0, ms = 5, color = "C0", label = rdspl(abs(d)),)

ax.legend()
ax.set_ylabel("log(⟨t⟩)")
ax.set_xlabel("log(ε)")
ax.set_title("Towel map")
fig.tight_layout()
