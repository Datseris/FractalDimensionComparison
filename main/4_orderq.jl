# %% Sensititivy to trajectory length
using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

# %%
N = 1*10^5
systems = [:koch, :henon_chaotic]
slabels = ["Koch", "Hénon"]
qs = 2:4
Cmethod = "standard" # bueno or standard. Decides εmax for correlation sum.

eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for data in systems
    for q in qs
        qH = qC = q
        @show (q, data)
        # Here we simply pack all parameters into a dictionary
        # (other parameters are (probably) globals)
        params = @strdict N qH qC data
        if data == :standardmap_chaotic
            params["k"] = 1.0
        # elseif data == :henon_chaotic
            # params["z"] = -4
        end
        params["theiler"] = 0
        if Cmethod ≠ "standard"
            params["Cmethod"] = Cmethod
        end

        # This is the main call that calculates everything
        output, s = produce_or_load(
            datadir("main"), params, make_C_H;
            prefix = string(data), suffix = "jld2", force = false,
            ignores = ["data"], storepatch = false
        )

        @unpack eH, eC, H, C = output
        push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
    end
end

legendtitle = "impact of order \$q\$"
labels = [s*" \$q=$(q)\$" for s in slabels for q in qs]

fig, axs = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle; 
    qH = "q", qC = "q", tol = 0.25, 
    offsets = range(0; length = 6, step = 1.5),
    dimension_fit_C = FractalDimension.linear_regression_fit_glm,
)

wsave(plotsdir("paper", "orderq"), fig)
