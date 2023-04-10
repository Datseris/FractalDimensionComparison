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
