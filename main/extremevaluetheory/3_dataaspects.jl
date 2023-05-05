using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

p = 0.99

# 1st subplot: length
Ns = reverse([500, 1000, 5000, 10000, 50000, 100_000])
labels_length = ["N=$(N)" for N in Ns]
data = :lorenz96_chaotic

Dlocs_length = Vector{Float64}[]
for N in Ns
    params = @strdict data p N
    if data == :lorenz96_chaotic
        params["D"] = 8
    end
    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs_length, Δloc)
end

# 2nd subplot: dimension
N = Int(1e5)
system = :henonheiles_chaotic # change it to whatever system you want
data = :embed_system
embedding_ds = 3:8
labels_embed = ["d=$(d)" for d in embedding_ds]
Dlocs_embed = Vector{Float64}[]
for d in embedding_ds
    params = @strdict data p N d system
    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs_embed, Δloc)
end

# 3rd subplot: sampling time
data = :roessler_chaotic # change it to whatever system you want
times = [0.001, 0.1, 10.0]
labels_times = ["Δt=$(t)" for t in times]
N = Int(1e5)
Dlocs_times = Vector{Float64}[]

for Δt in times
    params = @strdict data p N Δt
    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs_times, Δloc)
end

# %% plotty plot
# Special figure
fig = Figure(resolution = (figwidth, 2figheight))
axs = [Axis(fig[i, 1]) for i in (2, 4, 6)]
legend_pos = [fig[i, 1] for i in (1, 3, 5)]

evtplot!(axs[1], Dlocs_length, labels_length, "varying length";
    legend_position = legend_pos[1],
    upperlim = 16, lowerlim = 2,
    expected = fill(6.91, 6), cutoffs = fill(8, 6),
)

evtplot!(axs[2], Dlocs_embed, labels_embed, "input dimensionality (delay embedding)";
    legend_position = legend_pos[2],
    upperlim = 6, lowerlim = 1.5,
    expected = fill(3, 6), cutoffs = fill(3, 6),
)

evtplot!(axs[3], Dlocs_times, labels_times, "sampling time";
    legend_position = legend_pos[3],
    upperlim = 4, lowerlim = 1,
    xpositions = [1, 3, 5],
    expected = fill(1.9, 6), cutoffs = fill(3, 6),
)

for i in 1:2
    hidexdecorations!(axs[i]; grid = false)
end

fig

wsave(plotsdir("paper", "evt_dataaspects"), fig)
