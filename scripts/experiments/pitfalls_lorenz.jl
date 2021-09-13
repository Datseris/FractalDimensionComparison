using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))

# %%
tr = Data.lorenz96_nonstationary(;N=1e4, D = 4)
x,y,z = columns(tr)
figure()
plot3D(x, y, z; lw = 1.0, color = "C0", alpha = 0.2)
scatter3D(x, y, z; lw = 1.0, color = "C0", alpha = 0.2)
# compare this with actual attractor
tr = Data.lorenz96_chaotic(;N=1e4, D = 4)
x,y,z = columns(tr)
figure()
plot3D(x, y, z; lw = 1.0, color = "C0", alpha = 0.2)
scatter3D(x, y, z; color = "C0", alpha = 0.2)

# %%
datas = Symbol[]
labels = String[]

N = 1e5

push!(datas, :lorenz96_nonstationary)
push!(labels, "non-stationary Lorenz96")

push!(datas, :lorenz96_chaotic)
push!(labels, "standard Lorenz96")


qH = 1
qC = 2
# Calculate values for H, C
es, Hs, Cs = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
for i in 1:length(datas)
    data = datas[i]
    params = @strdict N qH qC data
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        prefix = string(data), suffix = "jld2", force = false,
        ignores = ["data"],
    )
    @unpack e, H, C = output
    push!(es, e); push!(Hs, H); push!(Cs, C)
end

# Do the actual plot
legendtitle = "inappropriate data"
fig, axs = mainplot(Hs, Cs, es, labels, legendtitle; qH, qC, tol = 0.25)
