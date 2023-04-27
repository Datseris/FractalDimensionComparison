# %% Sensititivy to trajectory length
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

N = Int(1e5)
Ds = range(4; step = 2, length = 6)

qH = 1
qC = 2
data = :lorenz96_chaotic

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(Ds)
    # These are the parameters that change between the different elements of the plot
    D = Ds[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict N D qH qC data

    # This is the main call that calculates everything
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = true,
        ignores = ["data"],
    )
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end

# Do the actual plot
labels = ["\$D=$(D)\$" for D in Ds]
legendtitle = "impact of dimensionality \$D\$ (Lorenz96 system)"

fig, axs = mainplot(
    Hs, Cs, es, labels, legendtitle;
    qH, qC, tol = 0.25,
    # dimension_fit = FractalDimension.logarithmic_corrected_fit_lsqfit
)

# wsave(plotsdir("paper", "dimensionality"), fig)
