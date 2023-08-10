# %% Sensititivy to trajectory length
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

N = Int(10^5)
systems = [:koch, :henon_chaotic]
slabels = ["Koch", "Hénon"]
qs = 2:4
# qs = Real[2, 2.01, 4] # test that q=2 and q=2.01 give identical results
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
        elseif data == :lorenz96_chaotic
            params["D"] = 8
        end
        params["theiler"] = 0
        if Cmethod ≠ "standard"
            params["Cmethod"] = Cmethod
        end

        # This is the main call that calculates everything
        output = produce_or_load_C_H(params, data; force = false)

        @unpack eH, eC, H, C = output
        push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)
    end
end

legendtitle = "impact of order q"
labels = [s*" q=$(q)" for s in slabels for q in qs]

fig = mainplot(
    Hs, Cs, eHs, eCs, labels, legendtitle;
    qH = "q", qC = "q", tol = 0.25,
    offsets = range(0; length = 6, step = 1.5),
    dimension_fit_C = linear_regression_fit_linalg,
    region_choice = :last,
)

wsave(plotsdir("paper", "orderq"), fig)

# %% Analysis of where the correlation sum slopes change and how much are they

fig = Figure()

for j in (5, 6) # last two is the two slopes
    C = Cs[j]
    eC = eCs[j]
    ax = Axis(fig[j-4, 1])

    i = findfirst(c -> c > 0, C)
    if !isnothing(i)
        x, y = log.(eC)[i:end], log.(C)[i:end]
    else
        error()
    end

    lrs, tangents = linear_regions(x, y)

    for r in lrs
        scatterlines!(ax, x[r], y[r])
    end

    println("Data $j")
    for k in eachindex(lrs)
        @show (lrs[k], tangents[k])
    end
end

fig