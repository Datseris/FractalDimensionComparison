# %% Demonstration of pitfalls/invalid sets
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

datas = Symbol[]
labels = String[]

N = Int(1e4)
p = 0.98
qH = 2
qC = 2
z = 1

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

Dlocs = Vector{Float64}[]
for i in eachindex(datas)
    # These are the parameters that change between the different elements of the plot
    data = datas[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    if string(data)[1:3] != "exp"
        params = @strdict N data
    else
        name = labels[i]
        params = @strdict data name
    end
    if data == :lorenz96_chaotic
        D = 32
        N = Int(1e5)
        @pack! params = D, N
    end
    if data == :ksiva
        N = Int(1e5)
        @pack! params = N
    end

    params_evt = copy(params)
    params_c = copy(params)
    @pack! params_evt = p
    @pack! params_c = qH, qC, z

    # This is the main call that calculates everything
    output = produce_or_load_EVT(params_evt, data; force = false)
    @unpack Δloc = output
    if data == :lorenz96_chaotic || data == :ksiva
        Δloc ./= 4
    end

    push!(Dlocs, Δloc)
end

# Do the actual plot
legendtitle = "extreme cases, p = $(p)"
fig = evtplot(Dlocs, labels, legendtitle;
    cutoffs = [6, 3, Inf, Inf, Inf, Inf],
    upperlim = 8, lowerlim = 0.5,
)

display(fig)
# wsave(plotsdir("paper", "evt_extreme"), fig)
