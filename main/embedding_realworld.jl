using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
N = Int(1e5)
data = :embed_system
system = :experimental_data
name = "vostok"
legendtitle = "delay embedding real world data: $(name)"
embedding_ds = 3:8

# params for H,C
qH = 2
qC = 2

# params for EVT
p = 0.95
estimator = :exp

eHs, eCs, Hs, Cs, Dlocs = [Vector{Float64}[] for i in 1:5]

for d in embedding_ds
    params = @strdict N d data system name
    params_evt = copy(params)
    @pack! params_evt = p, estimator
    params_hc = copy(params)
    @pack! params_hc = qH, qC

    output = produce_or_load_C_H(params_hc, data; force = false)
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)

    output = produce_or_load_EVT(params_evt, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)
end

labels = ["d=$(d)" for d in embedding_ds]

fig = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tolH = 0.25, tolC = 0.25,
    offsets = reverse(range(0; length = 6, step = 1.5)),
    region_choice = :last,
    dimension_fit_C = linear_regression_fit_linalg,
    resolution = (figwidth, 1.5figheight),
)

axevt = Axis(fig[3,1])
evtplot!(axevt, Dlocs, labels, legendtitle;
    upperlim = 5, lowerlim = 1.5,
    legend_position = nothing,
    # cutoffs = fill(3, 6),
    # expected = fill(1.9, 6),
)
display(fig)

wsave(plotsdir("paper", "embedding_experimental_$(name)"), fig)
