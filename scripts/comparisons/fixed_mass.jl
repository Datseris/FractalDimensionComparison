using DrWatson
@quickactivate :FractalDimensionComparison

# Input data
datasets = Vector(undef, 6)
labels = Vector{String}(undef, 6)

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

# Meta-parameters
N = Int(1e5)
qH = qC = 2
max_j = 20
start_j = 3

# Compute and store quantities
eCs, Cs, rs, ys = [Vector{Float64}[] for i in 1:4]
for i in eachindex(datasets)
    data = datasets[i]
    params = @strdict N data
    if data == :lorenz96_chaotic
        params["D"] = D
    end
    # Corr sum
    paramsC = copy(params)
    @pack! paramsC = qH, qC
    output = produce_or_load_C_H(paramsC, data; force = false)
    @unpack eC, C = output
    push!(eCs, eC); push!(Cs, C)
    # Fixed mass (using `DrWatson.produce_or_load`)
    paramsF = copy(params)
    @pack! paramsF = max_j, start_j
    output, = produce_or_load(paramsF, datadir("fixedmass");
            filename = params -> savename(params; ignores = ["data"]),
            prefix = string(data)*"_fixedmass", suffix = "jld2", storepatch = false,
            force = false,
        ) do config
        @unpack max_j, start_j, data = paramsF
        @info "Producing data..."
        data_producing_function = getfield(Data, data)
        @time X = data_producing_function(; dict2ntuple(params)...)
        @info "Calculating fixed mas..."
        r, y = fixedmass_correlationsum(X, max_j; start_j)
        @pack! params = r, y
        return params
    end
    push!(rs, output["r"]); push!(ys, output["y"])
end

# %% Plot
using FractalDimensionComparison: rdspl
legendtitle = "corrsum vs. fixed mass"
# labels = ["N=$(N)" for N in Ns]
fig, axs = axesgrid(2, 1; sharex = false)
offsets = zeros(length(Cs))
llines = []
for j in eachindex(Cs)
    z = offsets[j]
    # corrsum
    C = Cs[j]
    eC = eCs[j]
    i = findfirst(c -> c > 0, C)
    isnothing(i) && error("All C is 0")
    x, y = log.(eC)[i:end], log.(C)[i:end]
    region, d = linear_region(x, y; warning = false)
    Δ, Δ05, Δ95 = linear_regression_fit_linalg(x[region], y[region])
    Clabel = "$(rdspl(Δ05))-$(rdspl(Δ95))"
    l = lines!(axs[1], x, y .+ z; alpha = 0.9)
    push!(llines, l)
    scatter!(axs[1], x[[region[1], region[end]]], y[[region[1], region[end]]] .+ z;
        label = Clabel, alpha = 0.75
    )
    # fixed mass
    x, y = rs[j], ys[j]
    region, d = linear_region(x, y; warning = false)
    Δ, Δ05, Δ95 = linear_regression_fit_linalg(x[region], y[region])
    Flabel = "$(rdspl(Δ05))-$(rdspl(Δ95))"
    lines!(axs[2], x, y .+ z; alpha = 0.9)
    scatter!(axs[2], x[[region[1], region[end]]], y[[region[1], region[end]]] .+ z;
        label = Flabel, alpha = 0.75
    )

end

axs[1].ylabel = L"\log(C_%$(qC))"
axs[1].xlabel = L"\log(\varepsilon)"
axs[2].xlabel = L"\langle \log \left( r(j)\right)\rangle,\quad j \in (%$start_j, %$max_j)"
axs[2].ylabel = L"\Psi(j)"

for ax in axs
    ax.xticks = WilkinsonTicks(8; k_min = 6)
    axislegend(ax; labelsize = 20, nbanks = 3, patchsize=(10f0,20), position = :lt, groupgap = 8)
end

# Make the informative legend
leg = Legend(fig[0, :], llines, labels, legendtitle;
    nbanks=3, patchsize=(40f0, 20),
)
rowgap!(fig.layout, 10)
space_out_legend!(leg)
space_out_legend!(leg) # not sure why I have to trigger this twice...

display(fig)

wsave(plotsdir("fixedmass", "different_systems_$(N)_max_j=$(max_j)"), fig)
wsave(plotsdir("paper", "fixedmass"), fig)
