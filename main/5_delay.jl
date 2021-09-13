using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))

N = Int(1e5)

qH = 1
qC = 2

system = :henonheiles_chaotic
data = :embed_system
embedding_ds = 2:7
Cmethod = "standard"

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for d in embedding_ds
    params = @strdict N d qH qC data system
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end

    params["z"] = -1

    # This is the main call that calculates everything
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = false,
        ignores = ["data"], tag = false,
    )
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end

# Do the actual plot
labels = ["\$d=$(d)\$" for d in embedding_ds]
legendtitle = "delay embedding (Hénon-Heiles, Δ ≈ 3)"

fig, axs = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle; 
    qH, qC, tol = 0.25, 
)

wsave(plotsdir("paper", "embedding_$(system)"), fig)
