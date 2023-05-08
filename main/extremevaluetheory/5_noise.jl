using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

N = Int(10^5)
p = 0.99

ds = range(2; step = 1, length = 6)

datas = fill(:roessler_chaotic, 3)
append!(datas, fill(:roessler_dynamic, 2))
push!(datas, :roessler_rounded)
# imo a 10% dynamic noise would correspond to a η≈2 for Rössler,
# as the regularisation happens afterwards
parameters = ["η" => x for x in [0, 0.05, 0.1, 1, 2, 0]]

Dlocs = Vector{Float64}[]

for i in eachindex(datas)
    data = datas[i]
    params = @strdict data p N
    push!(params, parameters[i])

    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)
end

labels = ["none", "5% additive", "10% additive", "5% dynamic", "10% dynamic", "rounding"]
legendtitle = "noise"

fig = evtplot(Dlocs, labels, legendtitle;
    upperlim = 4, lowerlim = 0.5, cutoffs = fill(3, 6),
    expected = fill(1.9, 6),
)

display(fig)

wsave(plotsdir("paper", "evt_noise"), fig)
