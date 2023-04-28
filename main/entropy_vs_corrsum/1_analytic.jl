# %% Comparison with analytically resolved models
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

datas = Vector(undef, 6)
labels = Vector{String}(undef, 6)

N = Int(1e5)

datas[1] = :roessler_periodic
labels[1] = "periodic"

# datas[2] = :torus2
# labels[2] = "2-torus"
datas[2] = :henonheiles_quasi
labels[2] = "quasiperiodic"

datas[3] = :koch
labels[3] = "Koch snowflake"

datas[4] = :kaplanyorke_map
labels[4] = "Kaplan-Yorke map"
# datas[4] = :brownian_motion
# labels[4] = "3D Brownian motion"

datas[5] = :uniform_sphere
labels[5] = "sphere"

datas[6] = :standardmap_chaotic
labels[6] = "SM (uniform)"

qH = 0
qC = 2
Cmethod = "standard" # bueno or standard. Decides εmax for correlation sum.

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(datas)
    # These are the parameters that change between the different elements of the plot
    data = datas[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict N qH qC data
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end
    if data == :koch
        params["maxk"] = 7
        delete!(params, "N")
    end
    if data == :brownian_motion
        params["theiler"] = 0
    end

    # This is the main call that calculates everything
    output = produce_or_load_C_H(params, data; force = true)
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end


# Do the actual plot
legendtitle = "analytically known Δ"

fig = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tol = 0.25,

    # For this plot we use standard regression fit because
    # the estimate is already so accurate, we don't want to
    # have the unecessary larger confidence intervals from the
    # logarithmic correction
    dimension_fit_C = linear_regression_fit_linalg,
)

wsave(plotsdir("paper", "analytic"), fig)
