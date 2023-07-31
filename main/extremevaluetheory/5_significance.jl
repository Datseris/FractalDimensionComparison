using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

using HypothesisTests, Distributions

# Arbitrarily vary the chosen systems or p.
# We do not save data with `produce_or_load` because they occupy several GBs...

N = Int(5e4)
ps = [0.95, 0.99, 0.999]
# ps = [0.9, 0.95, 0.99]

data = :experimental_data
label = "electrochemical 1"
name = "electrochemical1"
estimator = :exp
TestType = ExactOneSampleKSTest

function significance_plot!(fig, Es, nrmses, pvalues, sigmas, xis, p)
    # Choose indices by first sorting by pvalue
    sperm = sortperm(nrmses)
    js = (sperm[1], sperm[length(sperm)÷2], sperm[end])

    for i in 1:3
        ax = Axis(fig[i, 2])
        j = js[i] # rand(eachindex(Es))
        E = Es[j]
        s, x, pval = sigmas[j], xis[j], pvalues[j]
        gpd = GeneralizedPareto(0, s, x)
        nbins = 25

        hist!(ax, E; bins = nbins, normalization = :pdf)
        xrange = range(0, maximum(E); length = 100)
        nrmse = nrmses[j]
        t = "p-value = $(round(pval;digits = 3))\nnrmse = $(round(nrmse;digits = 3))"
        lines!(ax, xrange, pdf.(gpd, xrange); label = t, color = :black)
        Makie.text!(ax, 1, 1; text = t, align = (:right, :top), space = :relative, offset = (-4, -4))
        hideydecorations!(ax)
        ylims!(ax, 0, nothing)
        xlims!(ax, 0, nothing)
    end

    axh = Axis(fig[:, 1])

    bins = 0.0:0.05:1.0
    hist!(axh, pvalues; bins, label = "p-value, m = $(round(median(pvalues); digits = 2))", normalization = :pdf, color = COLORS[3])
    hist!(axh, clamp.(nrmses, 0, prevfloat(1.0)); bins, label = "NRMSE, m = $(round(median(nrmses); digits = 2))", normalization = :pdf, color = (COLORS[4], 0.8))
    axislegend("p = $(p)")
    xlims!(axh, 0, 1)
    ylims!(axh, 0, nothing)
    axh.xticks = ([0, 0.5, 1], ["0", "0.5", "1"])
    hideydecorations!(axh)
end

fig = Figure(resolution = (figwidth, 2figheight))

for (i, p) in enumerate(ps)
    gl = GridLayout(fig[i, :])

    params = @strdict data p N estimator TestType
    if data == :lorenz96_chaotic
        params["D"] = 8
    end
    if data == :experimental_data
        params["name"] = name
    end

    # params, = produce_or_load(params, datadir("evt");
    #     filename = params ->
    #         savename(params; ignores = ["data"], val_to_string = string, allowedtypes = (Real, Symbol, String, Type)),
    #     prefix = "evt_signif_"*string(data), suffix = "jld2", storepatch = false,
    #     force = false) do params

        data_producing_function = getfield(Data, data)
        X = data_producing_function(; dict2ntuple(params)...)
        Es, nrmses, pvalues, sigmas, xis = extremevaltheory_gpdfit_pvalues(X, p; estimator, TestType)

        @pack! params = Es, nrmses, pvalues, sigmas, xis
    #     return params
    # end

    @unpack Es, nrmses, pvalues, sigmas, xis = params

    significance_plot!(gl, Es, nrmses, pvalues, sigmas, xis, p)
end

figuretitle!(fig, "significance of EVT: $(label)")
display(fig)

wsave(plotsdir("paper", "evt_significance_$(name)_$(estimator)"), fig)
