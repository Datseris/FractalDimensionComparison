using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
N = Int(1e5)

system = :henonheiles_chaotic
legendtitle = "delay embedding + 5% noise"
data = :embed_system
embedding_ds = 3:8
η = 0.05 # noise level, percentage, from standardized timeseries
p = 0.95
estimator = :exp

Dlocs = []

for d in embedding_ds
    params = @strdict N d data system η p estimator


    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)
end


labels = ["d=$(d)" for d in embedding_ds]
legendtitle = "delay embedding + 5% noise"

fig = evtplot(Dlocs, labels, legendtitle;
    upperlim = 5, lowerlim = 1.5,
    # cutoffs = fill(3, 6),
    # expected = fill(1.9, 6),
)

# %% no noise but experimental

legendtitle = "delay embedding real world data"
data = :embed_system
system = :experimental_data
name = "nifty50"
embedding_ds = 3:8
p = 0.98
estimator = :exp

eHs, eCs, Hs, Cs, Dlocs = [Vector{Float64}[] for i in 1:5]

for d in embedding_ds
    params = @strdict N d data system p estimator name

    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)
end


labels = ["d=$(d)" for d in embedding_ds]

fig = evtplot(Dlocs, labels, legendtitle;
    upperlim = 5, lowerlim = 1.5,
    # cutoffs = fill(3, 6),
    # expected = fill(1.9, 6),
)
