# %% Demonstration of pitfalls/invalid sets
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

datas = Symbol[]
labels = String[]

N = Int(1e4)
dt = 1.0

push!(datas, :lorenz96_nonstationary)
push!(labels, "non-stat. Lorenz96")

push!(datas, :roessler_sphere)
push!(labels, "Periodic+3Dsphere")

push!(datas, :experimental_data)
push!(labels, "vostok")

push!(datas, :experimental_data)
push!(labels, "nifty50")

push!(datas, :lorenz96_chaotic)
push!(labels, "Lorenz96 (D=32)")

push!(datas, :ksiva)
push!(labels, "K-Sivashinsky")

qH = 2
qC = 2
Cmethod = "standard" # bueno or standard. Decides εmax for correlation sum.

# Calculate values for H, C
eHs, eCs, Hs, Cs = [Vector{Float64}[] for i in 1:4]
for i in eachindex(datas)
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
    params["z"] = 0
    if data == :lorenz96_chaotic
        N = Int(1e5)
        D = 32
        prism = 32
        use_boxed = false
        params["w"] = 0
        e = MathConstants.e .^ (-1:0.2:3)
        @pack! params = D, N, prism, use_boxed, e
    end
    if data == :ksiva
        N = Int(1e5)
        prism = 100
        use_boxed = false
        params["w"] = 0
        e = MathConstants.e .^ (-1:0.1:3)
        compute_H = false
        @pack! params = N, prism, use_boxed, e, compute_H
    end

    # This is the main call that calculates everything
    output = produce_or_load_C_H(params, data; force = i > 4)

    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
end


# Do the actual plot
legendtitle = "extreme cases"

fig = mainplot(Hs, Cs, eHs, eCs, labels, legendtitle; qH, qC, tol = 0.25,
    dimension_fit_C = linear_regression_fit_linalg,
)

display(fig)
wsave(plotsdir("paper", "extreme"), fig)
