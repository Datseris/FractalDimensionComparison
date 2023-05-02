using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
N = Int(1e5)

qH = 2
qC = 2

system = :towel_chaotic
legendtitle = "delay embedding (Towel map, Δ ≈ 2)"
data = :embed_system
embedding_ds = 3:8
Cmethod = "standard"
theiler = τ = 1 # use nothing to estimate via Mutual Information first minimum

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for d in embedding_ds
    params = @strdict N d qH qC data system
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end

    params["z"] = -1
    params["theiler"] = params["τ"] = theiler

    # This is the main call that calculates everything
    output = produce_or_load_C_H(params, data; force = false)

    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end

# Do the actual plot
labels = ["d=$(d)" for d in embedding_ds]

fig = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tol = 0.25,
    offsets = reverse(range(0; length = 6, step = 1.5)),
)

wsave(plotsdir("paper", "embedding_$(system)"), fig)
