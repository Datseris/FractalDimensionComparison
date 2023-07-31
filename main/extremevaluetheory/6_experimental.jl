using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

names = String[
    "electrochemical1", "electrochemical2", "shinriki1",
    "roessler_embed",  "doublependulum", "ceps",
]
labels = String[
    "electroch. 1", "electroch. 2", "Shinriki",
    "Rössler Net",  "Double pend.", "ECG IBI",
]

estimator = :exp
p = 0.98
Dlocs = Vector{Float64}[]

for i in 1:length(labels)
    # These are the parameters that change between the different elements of the plot
    data = :experimental_data
    name = names[i]
    params = @strdict p data name estimator

    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)
end

legendtitle = "experimental data, p = $(p)"
fig = evtplot(Dlocs, labels, legendtitle;
    upperlim = 5, lowerlim = 0.5,
)

display(fig)
wsave(plotsdir("paper", "evt_experimental_$(estimator)"), fig)
