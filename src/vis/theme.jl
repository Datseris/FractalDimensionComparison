using CairoMakie
import Downloads

path = joinpath(@__DIR__, "plottheme.jl")

try
    Downloads.download(
        "https://raw.githubusercontent.com/Datseris/plottheme/main/plottheme.jl",
        path
    )
catch
end
if isfile(path)
    include(path)
end

# Specific to this paper:
include("mainplot.jl")

# Figure sizes corresponding to default (1 column) figures
# with three rows (one is legend)
figwidth = 800
figheight = 600

# Because we use DrWatson's `@quickactivate :MyProjectName`
function __init__()
    set_theme!(default_theme)
    update_theme!(;
        resolution = (figwidth, figheight),
        Lines = (cycle = Cycle([:color, :linestyle], covary = true), ),
        Scatter = (cycle = Cycle([:color, :marker], covary = true), markersize = 20),
    )
end

# Extend DrWatson's `wsave`
DrWatson._wsave(s, fig::Makie.Figure) = Makie.save(s*".png", fig; px_per_unit = 2)
