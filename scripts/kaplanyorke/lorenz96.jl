#=
This file uses the Kaplan Yorke dimension to provide high accuracy estimates for
the fractal dimension for the Lorenz-96 system.
=#
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
using DynamicalSystems, Statistics, OrdinaryDiffEq

Dmax = 13
ics = 10 # how many random initial conditions to obtain
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9, maxiters = typemax(Int))
Ds = 5:Dmax
K = [Float64[] for _ ∈ 1:length(Ds)]
Λ = [Vector{Float64}[] for _ ∈ 1:length(Ds)]
T = 100000.0
F = 8.0

for (i, D) ∈ enumerate(Ds)
    @show D
    for j ∈ 1:ics
        ds = Systems.lorenz96(D; F = F)
        @time λs = lyapunovs(ds, T; Ttr = 100.0, diffeq...)
        push!(Λ[i], λs)
        push!(K[i], kaplanyorke_dim(λs))
    end
end

name = savename("lorenz96", @dict(Dmax, ics, T, F), "jld2")
result = @strdict Ds ics K Λ
@tagsave(datadir(name), result)
