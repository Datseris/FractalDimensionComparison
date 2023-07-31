using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

N = Int(1e5)
p = 0.99
estimator = :exp

datas = Vector(undef, 6)
labels = Vector{String}(undef, 6)

datas[1] = :henonheiles_quasi
labels[1] = "quasiperiodic"

datas[2] = :koch
labels[2] = "Koch snowflake"

datas[3] = :lorenz96_chaotic
labels[3] = "1/4 Lorenz96"

datas[4] = :henon_chaotic
labels[4] = "Hénon map"

datas[5] = :standardmap_chaotic
labels[5] = "standard map"

datas[6] = :coupled_logistics
labels[6] = "1/4 Coupled logistics"

Dlocs = Vector{Float64}[]
Clocs = Vector{Float64}[]

for i in 1:length(datas)
    # These are the parameters that change between the different elements of the plot
    data = datas[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict data N

    if data == :koch
        params["maxk"] = 7
        delete!(params, "N")
    elseif data == :lorenz96_chaotic
        params["D"] = 8
    end

    # This is the main call that calculates everything
    paramsevt = deepcopy(params)
    paramsevt["estimator"] = estimator
    paramsevt["p"] = p
    output = produce_or_load_EVT(paramsevt, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)

    output = produce_or_load_pointwise(params, data; force = false)
    Cloc = output["Ds"]
    Csums = output["Cs"]
    L = length(Csums[1])÷2
    first_entries = [C[1] for C in Csums]
    last_entries = [C[end] for C in Csums]
    mid_entries = [C[L] for C in Csums]
    @show i, median(mid_entries)

    @show count(iszero, Cloc)
    filter!(!iszero, Cloc)
    push!(Clocs, Cloc)

end

for j in (3, 6)
    Dlocs[j] = Dlocs[j] ./ 4
    Clocs[j] = Clocs[j] ./ 4
end


fig = evtplot(Dlocs, labels, "comparison with pointwise";
    upperlim = 3, lowerlim = 0.5,
    side = :left, inner_legend = false,
    # cutoffs = [2, 3, 2, 2, 3, 2],
    # expected = [1, 2, 6.91/4, 1.26, 2.24, 8/4]
)

# Add the pointwise distributions
ax = content(fig[1,1])
for (i, Dloc) in enumerate(Clocs)
    vl = violin!(ax, fill(i, length(Dloc)), Dloc;
        color = (COLORS[i], 0.1), strokewidth = 2, strokecolor = COLORS[i], side = :right
    )
    horizontal = [i, i + 0.5 - 0.2/2]
    m = mean(Dloc)
    lines!(ax, horizontal, [m, m]; color = COLORS[i], linestyle = :dot)
end

display(fig)

wsave(plotsdir("paper", "evt_pointwise_$(estimator)"), fig)
