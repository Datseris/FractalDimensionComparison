using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using ComplexityMeasures # for histogram

data = :henon_chaotic
data_producing_function = getfield(Data, data)
X = data_producing_function(; N = 100_000)

xlims = (1.405, 1.41)
ylims = (-0.48, -0.33)
lims = (xlims, ylims)
Xsmall_idx = map(p -> xlims[1] < p[1] < xlims[2] && ylims[1] < p[2] < ylims[2], X)
Xsmall = X[Xsmall_idx]

# %%

fig = Figure(resolution = (1000, 600),; backgroundcolor = :transparent)
ax = Axis(fig[1,1], backgroundcolor = :transparent)

ms = Observable(2)
scatter!(ax, X[:, 1], X[:, 2]; markersize = ms, color = :black)
hidedecorations!(ax)
hidespines!(ax)
fig

wsave(plotsdir("summary", "henon_attractor.png"), fig)

# %%
# zoom-in
ax.limits = ((1.405, 1.41), (-0.5, -0.3))
fig = Figure(resolution = (600, 600),; backgroundcolor = :transparent)
ax = Axis(fig[1,1], backgroundcolor = :transparent)
scatter!(ax, Xsmall[:, 1], Xsmall[:, 2]; markersize = 6, color = :black)
hidedecorations!(ax)
hidespines!(ax)
fig

wsave(plotsdir("summary", "henon_attractor_zoom.png"), fig)


# %% entropy boxing figure
fig = Figure(resolution = (1000, 600),; backgroundcolor = :transparent)
ax = Axis(fig[1,1], backgroundcolor = :transparent)
scatter!(ax, Xsmall[:, 1], Xsmall[:, 2]; markersize = 6, color = :black)
ax.limits = lims

resize!(fig, (600, 600))

esx = [0.02, 0.005] ./ 12
esy = [0.04, 0.008]

for (i, (ex, ey)) ∈ enumerate(reverse(collect(zip(esx, esy))))
    p, b = probabilities_and_outcomes(ValueHistogram([ex, ey]), Xsmall)
    x = [a[1] for a in b]
    y = [a[2] for a in b]
    rects = [Rect(x[i], y[i], ex, ey) for i in eachindex(x)]
    vlines!(ax, minimum(x):ex:(maximum(x)+2ex); color = COLORS[i+2], linewidth = 0.75*i)
    hlines!(ax, minimum(y):ey:(maximum(y)+2ey); color = COLORS[i+2], linewidth = 0.75*i)

    poly!(ax, rects; color = (COLORS[i+2], 0.2))
end

hidedecorations!(ax)
hidespines!(ax)


fig

wsave(plotsdir("summary", "henon_attractor_boxed.png"), fig)


# %% histogram

# generate bins for histogram

ε = 0.01
p, b = probabilities_and_outcomes(ValueHistogram(ε), X)
x = [a[1] for a in b]
y = [a[2] for a in b]
p = p ./ maximum(p)

fig = Figure(resolution = (1000, 1000),; backgroundcolor = :transparent)
ax = Axis3(fig[1,1], backgroundcolor = :transparent)

# ax.limits=Rect3(Vec3f(0), Vec3f(1))
meshscatter!(ax, vec(Point3f.(x, y, 0.0)),
markersize=Vec3f.(2ε, 2ε, p), marker=Rect3(Vec3f(0), Vec3f(1)),
color = clamp.(p, 0, 0.25), alpha = 0.5
)

# %%
ax.azimuth = -π/3 - 0.1 #2.175530633326986
ax.elevation = 0.6
zlims!(ax, 0, 1)
hidedecorations!(ax)

fig

wsave(plotsdir("summary", "henon_attractor_density.png"), fig)


# %% Joint plot for top right corner
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using ComplexityMeasures # for histogram

N = Int(1e5)

datas = Vector(undef, 6)
labels = Vector{String}(undef, 6)

datas[1] = :koch
labels[1] = "Koch snowflake"

datas[2] = :henonheiles_quasi
labels[2] = "quasiperiodic"

datas[3] = :lorenz96_chaotic
labels[3] = "Lorenz96 (D=8)"

datas[4] = :henon_chaotic
labels[4] = "Hénon map"

datas[5] = :towel_chaotic
labels[5] = "Towel map"

datas[6] = :henonheiles_chaotic
labels[6] = "Hénon-Heiles"

# params for H,C
qH = 2
qC = 2
# params for EVT
p = 0.99
estimator = :exp

eHs, eCs, Hs, Cs, Dlocs = [Vector{Float64}[] for i in 1:5]

for data in datas
    params = @strdict N data
    params_evt = copy(params)
    @pack! params_evt = p, estimator
    params_hc = copy(params)
    @pack! params_hc = qH, qC

    output = produce_or_load_C_H(params_hc, data; force = false)
    @unpack eH, eC, H, C = output
    push!(eHs, eH); push!(Hs, H); push!(eCs, eC); push!(Cs, C)

    output = produce_or_load_EVT(params_evt, data; force = false)
    @unpack Δloc = output

    if data == :lorenz96_chaotic
        Δloc = Δloc/4
    end
    push!(Dlocs, Δloc)
end

# %%
fig = Figure(resolution = (1200, 850); backgroundcolor = :transparent)
cgl = GridLayout(fig[1, 1])
egl = GridLayout(fig[1, 2])

axh = Axis(cgl[1,1]; backgroundcolor = :transparent)
axc = Axis(cgl[2,1]; backgroundcolor = :transparent)
Makie.linkxaxes!(axh, axc)
Makie.hidexdecorations!(axh; grid = false)

legendtitle = "fractal dimensions of exemplary sets"

mainplot!(
    [axh, axc], Hs, Cs, eHs, eCs, labels, legendtitle;
    qH, qC, tolH = 0.25, tolC = 0.25,
    offsets = reverse(range(0; length = 6, step = 1.5)),
    region_choice = :last,
    dimension_fit_C = linear_regression_fit_linalg,
    legendpos = fig[0, :]
)

axevt = Axis(egl[1,1]; backgroundcolor = :transparent)
evtplot!(axevt, Dlocs, labels, legendtitle;
    upperlim = 4, lowerlim = 0.5,
    legend_position = nothing,
)
display(fig)

wsave(plotsdir("summary", "overview.png"), fig)
