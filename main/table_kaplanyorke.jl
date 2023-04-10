# %% Calculating Kaplan Yorke dimensions for as many systems as possible
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

# Time to evolve for reaching convergence
N = Int(1e6)
Ds = range(4; step = 2, length = 6)

function Δλ(params)
    @unpack N, D = params
    ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = 24.0)
    λs = lyapunovspectrum(ds, N; Ttr = 100, FractalDimension.Data.diffeq...)
    Δ = kaplanyorke_dim(λs)
    return @strdict N D λs Δ
end

for D in Ds
    params = @strdict N D
    file, s = produce_or_load(
        datadir("kaplanyorke"), params, Δλ;
        prefix = ("lorenz96_chaotic"), suffix = "jld2", force = true
    )
    println("For Lorenz-96:")
    println("D=$D, Δ=$(file["Δ"])")
end

# %% Show convergence explicitly
ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = 24.0)
tinteg = tangent_integrator(ds; FractalDimension.Data.diffeq...)
λs, t = ChaosTools.lyapunovspectrum_convergence(tinteg, N, 0.1, 100.0)
plot(t, λs)

# %% Rest of systems
systems = Symbol[
    :henon_chaotic,
    :towel_chaotic,
    :coupled_logistics, # D = 8
    :roessler_chaotic,
    :henonheiles_chaotic,
    # :lorenz96_chaotic, # D = 8, already calculated above
]
N = Int(1e6)

for s in systems # start producing
    params = @dict N

    ds = if s == :henon_chaotic
        Systems.henon([0.08152, 0.122717]; a = 1.4, b = 0.3)
    elseif s == :towel_chaotic
        Systems.towel([0.085, -0.121, 0.075])
    elseif s == :coupled_logistics
        params[:r] = 4.0
        u0 = range(0.1, 0.9; length = 8) .+ 0.0001
        ds = DiscreteDynamicalSystem(Data.logistics!, u0, [0.1, 3.0])
    elseif s == :roessler_chaotic
        ds = Systems.roessler([0.1, -0.2, 0.1]; c = 5.7, a = 0.2, b = 0.2)
    elseif s == :henonheiles_chaotic
        u0 = [0.0, -0.25, 0.42081, 0.0]
        ds = Systems.henonheiles(u0)
    end

    @show s
    function Δλ2(params)
        @unpack N = params
        λs = lyapunovspectrum(ds, N; Ttr = 100, FractalDimension.Data.diffeq...)
        @show λs
        Δ = kaplanyorke_dim(λs)
        return @strdict N λs Δ
    end

    file, s = produce_or_load(
        datadir("kaplanyorke"), params, Δλ2;
        prefix = string(s), suffix = "jld2", force = false
    )

    println("Δ=$(file["Δ"])")
end


# %% Try coupled tent maps
tr = Data.coupled_logistics(; D = 8, k = 0.1, r = 4.0)

using3D()
x,y,z= columns(tr)
scatter3D(x,y,z)
