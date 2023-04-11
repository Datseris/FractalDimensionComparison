# %% Fractal dimension of experimental datasets
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

names = String[
    "electrochemical1", "electrochemical2", "shinriki1",
    "roessler_embed",  "doublependulum",
    # "roessler_all", # doesn't make sense to use 28 dimensional data...
]
labels = String[
    "electroch. 1", "electroch. 2", "Shinriki",
    "Rössler Net",  "Double pend.",
    # "Rössler all",
]

all_ros_range = MathConstants.e .^ (-4:0.25:1)
Cmethod = "standard" # bueno or standard. Decides εmax for correlation sum.

qH = 1
qC = 2
# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(labels)
    # These are the parameters that change between the different elements of the plot
    data = :experimental_data
    name = names[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict qH qC data name
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end
    params["z"] = -1

    # This is the main call that calculates everything
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = false,
        ignores = ["data"], storepatch = false,
    )
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end


# Do the actual plot
legendtitle = "experimental data"

fig, axs = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tol = 0.25,
    offsets = range(0; length = 6, step = 2),
    # dimension_fit_C = FractalDimension.logarithmic_corrected_fit_lsqfit,
    dimension_fit_C = FractalDimension.linear_regression_fit_linalg,
)
wsave(plotsdir("paper", "experimental"), fig)
