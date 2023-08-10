using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

data1 = :lorenz96_chaotic
data2 = :towel_chaotic
data1_expected = 6.91/4
data2_expected = 2.2
data1_cutoff = 2
data2_cutoff = 3
estimator = :exp

# 1st subplot: vary p
N = Int(1e5)
title_ps = "varying p: fixed N = $(N)"
ps = [0.95, 0.999, 0.9999]
labels_ps = ["$(system), p=$(p)" for system in ("Lorenz96", "towel") for p in ps]
Dlocs_ps = Vector{Float64}[]

for data in (data1, data2)
    for p in ps
        params = @strdict data p N estimator
        if data == :lorenz96_chaotic
            params["D"] = 8
        end
        output = produce_or_load_EVT(params, data; force = false)
        @unpack Δloc = output
        if data ==:lorenz96_chaotic
            Δloc = Δloc ./ 4
        end
        push!(Dlocs_ps, Δloc)
    end
end

# 2nd subplot: fixed N(1-p)
set = 1000  # fixed product N p
ps = [0.99, 0.95, 0.9] # values of p to iterate over
labels_set = ["$(system), p=$(p)" for system in ("Lorenz96", "towel") for p in ps]
title_set = "varying p: fixed N(1-p) = $(round(Int, set))"
Dlocs_set = Vector{Float64}[]

for data in (data1, data2)
    for p in ps
        N = round(Int, set/(1 - p))
        @show N
        params = @strdict data p N estimator
        if data == :lorenz96_chaotic
            params["D"] = 8
        end
        output = produce_or_load_EVT(params, data; force = false)
        @unpack Δloc = output
        if data ==:lorenz96_chaotic
            Δloc = Δloc ./ 4
        end
        push!(Dlocs_set, Δloc)
    end
end


# %% plotty plot
# Special figure
fig = Figure(resolution = (figwidth, 1.25figheight))
axs = [Axis(fig[i, 1]) for i in (2, 4)]
legend_pos = [fig[i, 1] for i in (1, 3)]

evtplot!(axs[1], Dlocs_ps, labels_ps, title_ps;
    legend_position = legend_pos[1],
    upperlim = 4, lowerlim = 1,
    expected = [fill(data1_expected, 3)..., fill(2.2, 3)...],
    cutoffs = [fill(data1_cutoff, 3)..., fill(data2_cutoff, 3)...],
)

evtplot!(axs[2], Dlocs_set, labels_set, title_set;
    legend_position = legend_pos[2],
    upperlim = 4, lowerlim = 1,
    expected = [fill(data1_expected, 3)..., fill(data2_expected, 3)...],
    cutoffs = [fill(data1_cutoff, 3)..., fill(data2_cutoff, 3)...],
)

for i in 1:2
    hidexdecorations!(axs[i]; grid = false)
end

display(fig)

wsave(plotsdir("paper", "evt_quantile_$(estimator)"), fig)
