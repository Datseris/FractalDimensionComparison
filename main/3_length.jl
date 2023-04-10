# %% Sensititivy to trajectory length
using DrWatson
@quickactivate "FractalDimension"
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

D = 8
Ns = reverse([500, 1000, 5000, 10000, 20000, 50000])
labels = ["\$N=$(N)\$" for N in Ns]

qH = 1
qC = 2
Cmethod = "standard" # bueno or standard. Decides εmax for correlation sum.
# %%  Lorenz96

data = :lorenz96_chaotic

e = MathConstants.e .^ range(-2, 1; length = 16)

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(Ns)
    # These are the parameters that change between the different elements of the plot
    N = Ns[i]
    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict N D qH qC data
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end
    # For the length plot we use same range everywhere for consistency
    params["e"] = e
    # This is the main call that calculates everything
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = true,
        ignores = ["data"], tag = false, storepatch = false,
    )

    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end

# Do the actual plot
legendtitle = "impact of length \$N\$ (Lorenz96 system, \$D\$=$D)"

fig, axs = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tol = 0.25,
    # dimension_fit_C = FractalDimension.linear_regression_fit_glm,
    dimension_fit_C = FractalDimension.logarithmic_corrected_fit_lsqfit,
    offsets = range(0; length = 6, step = 1.5)
)
wsave(plotsdir("paper", "length_lorenz96"), fig)

# %% Repeat for towel map: 3D chaotic system
data = :towel_chaotic

eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(Ns)
    # These are the parameters that change between the different elements of the plot
    N = Ns[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict N qH qC data
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end

    # This is the main call that calculates everything
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = true,
        ignores = ["data"], storepatch = false,
    )
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end

# Do the actual plot
legendtitle = "impact of length \$N\$ (Towel map, \$D\$=3)"

fig, axs = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tol = 0.25,
    dimension_fit_C = FractalDimension.logarithmic_corrected_fit_lsqfit,
    offsets = 0:5
)
wsave(plotsdir("paper", "length_towel"), fig)

# %% Repeat for an experimental system
data = :experimental_data
name = "electrochemical2"
e = MathConstants.e .^ range(-3.5, 1; length = 16)

eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(labels)
    N = Ns[i]

    params = @strdict qH qC data name N
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end
    params["e"] = e

    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = true,
        ignores = ["data"], storepatch = false,
    )
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end

# Do the actual plot
legendtitle = "impact of length \$N\$ (experimental data)"

fig, axs = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tol = 0.25,
    dimension_fit_C = FractalDimension.logarithmic_corrected_fit_lsqfit,
    offsets = range(0; length = 6, step = 1.5)
)
wsave(plotsdir("paper", "length_experimental"), fig)
