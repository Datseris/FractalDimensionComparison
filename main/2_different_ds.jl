# %% Comparison with analytically resolved models
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

datasets = Vector(undef, 6)
labels = Vector{String}(undef, 6)

N = Int(1e5)

datasets[1] = :roessler_chaotic
labels[1] = "Rössler"

datasets[2] = :henonheiles_chaotic
labels[2] = "Hénon-Heiles"

datasets[3] = :lorenz96_chaotic
labels[3] = "Lorenz96 (D=8)"

datasets[4] = :henon_chaotic
labels[4] = "Hénon map"
datasets[5] = :towel_chaotic
labels[5] = "Towel map"
datasets[6] = :coupled_logistics
labels[6] = "Coupled logistic (D=8)"

qH = 1
qC = 2
Cmethod = "standard" # bueno or standard. Decides εmax for correlation sum.

eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(datasets)
    # These are the parameters that change between the different elements of the plot
    data = datasets[i]
    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict N qH qC data
    # For the two high dimensional systems we expand the ranges of
    # the ε-range, because by default it is too small
    if data == :lorenz96_chaotic
        params["D"] = 8
    end

    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end

    # This is the main call that calculates everything
    output = produce_or_load_C_H(params, data; force = false)
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end

legendtitle = "known dynamical rule"

fig = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tol = 0.25,
    offsets = [0, 0, -2, 0, 0, 0],
    # This chooses the kind of fit we do to the linear region.
    # The option `FractalDimension.logarithmic_corrected_fit_lsqfit`
    # is the method suggested by Sprott.
    # The alternative is `linear_regression_fit_linalg`,
    # which is standard linear regression.
    # By default, the correction is only done for the correlation sum method,
    # because this is what the Sprott paper says, but also because we found
    # (by experimentation) that this correction makes the entropy method
    # less accurate when used. For the correlation sum method the
    # accuracy of the corrected version is indeed better.

    # Feel free to experiment by altering any of these two:
    dimension_fit_C = logarithmic_corrected_fit_lsqfit,
    dimension_fit_H = linear_regression_fit_linalg,
)

wsave(plotsdir("paper", "different_systems_$N"), fig)
