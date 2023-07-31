using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
N = Int(1e5)

qH = 2
qC = 2

system = :henonheiles_chaotic
legendtitle = "delay embedding + 5% noise"
data = :embed_system
embedding_ds = 3:8
Cmethod = "standard"
η = 0.05 # noise level, percentage, from standardized timeseries

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for d in embedding_ds
    params = @strdict N d data system η
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end

    if system == :lorenz96_chaotic
        params["D"] = 8
    end
    # params["theiler"] = 1

    # This is the main call that calculates everything
    output = produce_or_load_C_H(params, data; force = false)

    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)

end

# Do the actual plot
labels = ["d=$(d)" for d in embedding_ds]

fig = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tolH = 0.25, tolC = 0.25,
    offsets = reverse(range(0; length = 6, step = 1.5)),
    region_choice = :last,
)

display(fig)

# wsave(plotsdir("paper", "embedding_$(system)"), fig)
