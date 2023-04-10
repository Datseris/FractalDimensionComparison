using PyPlot
include("style_colorscheme.jl")
using3D()

DrWatson._wsave(s, fig::Figure) = fig.savefig(s, dpi = 600, transparent = false)

"""
    rdspl(x, n = 2)
Round `x` to `n` sigdigits for display purposes.
"""
rdspl(x::Real, n = 2) = round(x, digits=n)
rdspl(x::AbstractVector, n = 2) = Tuple((round.(Float64.(x); sigdigits=n)))

PyPlot.rc("lines", lw = 3.0)
PyPlot.rc("errorbar", capsize = 6)
PyPlot.rc("axes", grid = true)
PyPlot.rc("grid", color = "0.75", alpha = 0.75)
PyPlot.rc("axes", axisbelow=true) # makes grid behind other plotting elements

PyPlot.rc("font", size = 24) # set default fontsize
PyPlot.rc("xtick", labelsize = 24)
PyPlot.rc("ytick", labelsize = 24)
PyPlot.rc("axes", labelsize = 28)
PyPlot.rc("legend";
    fontsize = 24, handletextpad = 0.5, handlelength = 2,
    labelspacing = 0.25, columnspacing = 1.0,
)
# PyPlot.rc("font", family = "Times New Roman") # Serif main font
PyPlot.rc("font", family = "DejaVu Sans") # sans main font
# PyPlot.rc("mathtext", rm = "sanserif", fontset="dejavusans") # sans math font
PyPlot.rc("mathtext", rm = "serif", fontset="cm") # serif math font

for z in ("x", "y")
    PyPlot.rc("$(z)tick.major", size = 7, width = 1.2)
    PyPlot.rc("$(z)tick.minor", size = 3, visible = false)
end

figx = 12 # default width correspoding to full text width
figy = 9  # default height corresponding to 1 row of plots
PyPlot.rc("figure", figsize = (figx, figy))
PyPlot.rc("savefig", dpi = 600, transparent = true, format = "png")

# set default color cycle
PyPlot.rc("axes", prop_cycle = matplotlib.cycler(color=COLORSCHEME))

if false # test color scheme
    fig = figure(figsize = (20, 15)) # show colors
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(223)
    ax5 = subplot(224)
    lw = 60
    L = length(COLORSCHEME)
    for (i, c) in enumerate(COLORS)
        chsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
        ax1.plot([0, 1], [0, 0] .+ i, color = c, lw = lw)
        ax1.set_title("color")
        ax2.plot([0, 1], [0, 0] .+ i, color = string(chsv[3]), lw = lw)
        ax2.set_title("brightness")
        ax3.plot([0, 1], [0, 0] .+ i, color = string(chsv[2]), lw = lw)
        ax3.set_title("saturation")
        x = 0:0.05:5π
        ax4.plot(x, cos.(x .+ i/2) .+ rand(length(x))/2; color=c, lw = 2)
        ax5.bar(collect(1:4) .+ (i-1)/L, 0.5rand(4) .+ 0.5, 1/L; color=c)
    end
    fig = tight_layout()
end

bbox = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 1.0)

"`add_identifiers!(fig = gcf(), axs = fig.get_axes(); xloc = 0.985, yloc = 0.975)`"
function add_identifiers!(fig = gcf(), axs = fig.get_axes(); xloc = 0.975, yloc = 1)
    bbox = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 1.0)
    for (i, ax) in enumerate(axs)
        l = ('a':'z')[i]
        try
            ax.text(xloc, yloc, "$(l)"; transform = ax.transAxes,
            bbox = bbox, zorder = 99, va = "top")
        catch e
            @warn string(e)
        end
    end
end
