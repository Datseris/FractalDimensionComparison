# %% Demonstration of pitfalls/invalid sets
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

datas = Symbol[]
labels = String[]

N = Int(1e4)
dt = 1.0

push!(datas, :lorenz96_nonstationary)
push!(labels, "non-stat. Lorenz96")

push!(datas, :roessler_sphere)
push!(labels, "Periodic+3Dsphere")

# push!(datas, :lorenz96_rounded)
# push!(labels, "rounded Lorenz96")

push!(datas, :experimental_data)
push!(labels, "vostok")

push!(datas, :experimental_data)
push!(labels, "nifty50")

qH = 1
qC = 2
Cmethod = "standard" # bueno or standard. Decides εmax for correlation sum.

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in 1:length(datas)
    # These are the parameters that change between the different elements of the plot
    data = datas[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    if string(datas[i])[1:3] != "exp"
        params = @strdict N qH qC data dt
    else
        name = labels[i]
        params = @strdict qH qC data name
    end
    if Cmethod ≠ "standard"
        params["Cmethod"] = Cmethod
    end
    params["z"] = 1

    # This is the main call that calculates everything
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = i == 3,
        ignores = ["data"], storepatch = false,
    )
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end


# Do the actual plot
legendtitle = "inappropriate data"

fig, axs = mainplot(Hs, Cs, eHs, eCs, labels, legendtitle; qH, qC, tol = 0.25)
wsave(plotsdir("paper", "inappropriate"), fig)
