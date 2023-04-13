#=
This module provides functions that generate the datasets
that are then used to calculate fractal dimensions.
The module is structured around functions which accept only keywords,
and always initialize and return a `StateSpaceSet`, by convention always given to
`standardize`. To add new datasets, simply create a new function here.
All functions use the `N` keyword, which is the length of the dataset,
and necessarily use keyword propagation (`kwargs...`) because of the
way they are called in `make_C_H.jl`.
Notice that you can easily also load data from CSV files or so, like
done in the function `experimental_data`. Also, embedding of the first timeseries
of any system can be done with the `embed_system` function.

Functions of this module are used like so:
```
data = :kaplanyorke_map # This must be a `Symbol`
data_producing_function = getfield(Data, data)
X = data_producing_function(; parameters...)
```
=#
module Data

using DrWatson
using Random, Statistics, DelimitedFiles, LinearAlgebra
using DynamicalSystemsBase, OrdinaryDiffEq, StochasticDiffEq

# default constants used if no other is given
const default_N = 10_000
const diffeq = (alg = Vern9(), reltol = 1e-12, abstol = 1e-12, maxiters = typemax(Int))
const diffeq_lowacc = (alg = Tsit5(), reltol = 1e-9, abstol = 1e-9, maxiters = typemax(Int))


# Kaplan Yorke
function kaplanyorke_map(; N = default_N, λ = 0.2, kwargs...)
    rng = Random.Xoshiro(1234)
    f(u, λ, t) = SVector((2u[1]) % 1.0 + 1e-15rand(rng), λ*u[2] + cospi(4u[1]))
    ds = DeterministicIteratedMap(f, SVector(0.15, 0.2), λ)
    tr, = trajectory(ds, N; Ttr = 1000)
    standardize(tr)
end

# Roessler
@inbounds function roessler_rule(u, p, t)
    du1 = -u[2]-u[3]
    du2 = u[1] + p[1]*u[2]
    du3 = p[2] + u[3]*(u[1] - p[3])
    return SVector{3}(du1, du2, du3)
end
function roessler_periodic(; N = default_N, Δt = 0.5, kwargs...)
    c = 3.0; a = 0.2; b = 0.2
    ds = CoupledODEs(roessler_rule, [0.1, -0.2, 0.1], [a, b, c]; diffeq)
    tr, = trajectory(ds, N*Δt; Δt, Ttr = 1000)
    standardize(tr)
end
function roessler_chaotic(; η=0, correlated=false, N = default_N, Δt = 0.5, kwargs...)
    c = 5.7; a = 0.2; b = 0.2
    ds = CoupledODEs(roessler_rule, [0.1, -0.2, 0.1], [a, b, c]; diffeq)
    tr, = trajectory(ds, N*Δt; Δt, Ttr = 1000)
    tr = standardize(tr)
    if η ≠ 0 # add static noise
        if !correlated
            seed = Random.MersenneTwister(5653435)
            for i in 1:length(tr)
                tr[i] += η*randn(seed, SVector{3})
            end
        else correlated
            a = arma_noise(N+1)
            for i in 1:length(tr)
                tr[i] = tr[i] .+ η .* a[i]
            end
        end
    end
    return tr
end
function roessler_embedding(; d = 3, kwargs...)
    x = roessler_chaotic(kwargs...)[:, 1]
    x = standardize(x)
    τ = estimate_delay(x, "mi_min")
    return embed(x, d, τ)
end
function roessler_rounded(; digits = 2, kwargs...)
    X = roessler_chaotic(; kwargs...)
    Z = [round.(x; digits = digits) for x in X]
    return StateSpaceSet(Z)
end
function roessler_dynoise(du,u, p, t)
    a, b, c, η = p
    du[1] = -u[2]-u[3]
    du[2] = u[1] + a*u[2]
    du[3] = b + u[3]*(u[1] - c)
end
# SDEProblem wants a specific function for the noise
# (This assumes diagonal noise, with the only noise acting on the
# 2nd variable)
function σ_roessler_dynoise(du, u, p, t)
    a, b, c, η = p
    du[1] = 0
    du[2] = η
    du[3] = 0
end
function roessler_dynamic(; N = default_N, Δt = 0.2, η=0.1, kwargs...)
    # using the SDE solve from DifferentialEquations here
    tspan = (0.0, N*Δt + 100)
    prob_sde_roessler = SDEProblem(
        roessler_dynoise, σ_roessler_dynoise, [0.1, -0.2, 0.1],
        tspan, [0.2, 0.2, 5.7, η]
    )
    saveat = 100.0:Δt:tspan[end]
    # I really don't have much experience with which solver to use
    sol = solve(prob_sde_roessler, SRA(); saveat, maxiters = typemax(Int))
    X = StateSpaceSet(hcat(sol.u...)')
    return standardize(X)
end
function roessler_lorenz(; N = default_N, kwargs...)
    tr1 = roessler_periodic(N = N÷2)
    tr2 = lorenz96_chaotic(; D = 4, N = N÷2)[:, 1:3]
    append!(standardize(tr1), standardize(tr2))
end
function roessler_sphere(; N = default_N, Δt = 0.5, kwargs...)
    ds = Systems.roessler([0.1, -0.2, 0.1]; c = 3.0, a = 0.2, b = 0.2; diffeq)
    tr1 = standardize(trajectory(ds, (N/2)*Δt; Δt, Ttr = 1000)[1])
    tr2 = standardize(uniform_sphere(; N = N÷2))
    return append!(tr1, tr2)
end
function roessler_nonstationary_eom(u, p, t)
    @inbounds begin
        a, b, dc = p
        du1 = -u[2]-u[3]
        du2 = u[1] + a*u[2]
        du3 = b + u[3]*(u[1] - u[4])
        du4 = dc
        return SVector{4}(du1, du2, du3, du4)
    end
end
function roessler_nonstationary(; N = default_N, Δt = 0.1, kwargs...)
    u0 = [-5.0, -5.0, 0.0, 3.0] # u0[4] is initial `c` value
    p0 = [0.2, 0.2, (20 - 3)/(N*Δt)]
    ds = ContinuousDynamicalSystem(roessler_nonstationary_eom, u0, p0; diffeq)
    tr, = trajectory(ds, N*Δt; Δt, Ttr = 0)
    return standardize(tr)
end

# Geometric sets
function torus2(; N = default_N, Δt = 0.1, ω = sqrt(3), R = 2.0, r = 1.0, kwargs...)
    function torus(u)
        θ, φ = u
        x = (R + r*cos(θ))*cos(φ)
        y = (R + r*cos(θ))*sin(φ)
        z = r*sin(θ)
        return SVector(x, y, z)
    end
    θs = range(0; step = 0.1, length = N)
    φs = ω .* θs
    return StateSpaceSet([torus(u) for u in zip(θs, φs)])
end
function koch(; maxk = 7, kwargs...)
    flakepoints = SVector{2}.([[0.0; 0.0], [0.5; sqrt(3)/2], [1; 0.0], [0.0; 0.0]])
    function innerkoch(points, maxk, α = sqrt(3)/2)
        Q = SMatrix(0, 1, -1, 0)
        for k = 1:maxk
            n = length(points)
            new_points = eltype(points)[]
            for i = 1:n-1
                p1, p2 = points[i], points[i+1]
                v = (p2 - p1) / 3
                q1 = p1 + v
                q2 = p1 + 1.5v + α * Q * v
                q3 = q1 + v
                push!(new_points, p1, q1, q2, q3)
            end
            push!(new_points, points[end])
            points = new_points
        end
        return points
    end
    kochpoints = innerkoch(flakepoints, maxk)
    kochdata = StateSpaceSet(unique!(kochpoints))
end
function uniform_sphere(; N = default_N, kwargs...)
    A = SVector{3, Float64}[]
    i = 0
    while i < N
        x = rand(SVector{3, Float64}) .- 0.5
        if norm(x) ≤ 1/2√3
            push!(A, x)
            i += 1
        end
    end
    return StateSpaceSet(A)
end

# Henon-Heiles
@inbounds function henonheiles_rule(u, p, t)
    du1 = u[3]
    du2 = u[4]
    du3 = -u[1] - 2.0*u[1]*u[2]
    du4 = -u[2] - (u[1]^2 - u[2]^2)
    return SVector(du1, du2, du3, du4)
end
function henonheiles_quasi(; N = default_N, Δt = 0.2, kwargs...)
    u0 = [0.0, 0.1, 0.5, 0.0]
    ds = CoupledODEs(henonheiles_rule, u0; diffeq)
    tr, = trajectory(ds, N*Δt; Δt, Ttr = 100)
    standardize(tr)
end
function henonheiles_chaotic(; N = default_N, Δt = 0.5, kwargs...)
    u0 = [0.0, -0.25, 0.42081, 0.0]
    ds = CoupledODEs(henonheiles_rule, u0; diffeq)
    tr, = trajectory(ds, N*Δt; Δt, Ttr = 100)
    standardize(tr)
end

# Standard map
@inbounds function standardmap_rule(x, par, n)
    theta = x[1]; p = x[2]
    p += par[1]*sin(theta)
    theta += p
    while theta >= 2π; theta -= 2π; end
    while theta < 0; theta += 2π; end
    while p >= 2π; p -= 2π; end
    while p < 0; p += 2π; end
    return SVector(theta, p)
end
function standardmap_chaotic(; N = default_N, k = 64.0, kwargs...)
    ds = DeterministicIteratedMap(standardmap_rule, [0.08152, 0.122717], [k])
    tr, = trajectory(ds, N; Ttr = 100)
    standardize(tr)
end

# Henon
henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
function henon_chaotic(; N = default_N, a = 1.4, b = 0.3, kwargs...)
    ds = DeterministicIteratedMap(henon_rule, [0.08152, 0.122717], [a, b])
    tr, = trajectory(ds, N; Ttr = 100)
    standardize(tr)
end

# Towel
function towel_rule(x, p, n)
    @inbounds x1, x2, x3 = x[1], x[2], x[3]
    SVector( 3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
    0.1*( (x2+0.35)*(1-2*x3) - 1 )*(1 - 1.9*x1),
    3.78*x3*(1-x3)+0.2*x2 )
end
function towel_chaotic(; N = default_N, η = 0, correlated = false, kwargs...)
    ds = DeterministicIteratedMap(towel_rule, [0.085, -0.121, 0.075])
    tr, = trajectory(ds, N; Ttr = 100)
    tr = standardize(tr)
    if η ≠ 0 # add static noise
        if !correlated
            seed = Random.MersenneTwister(5653435)
            for i in 1:length(tr)
                tr[i] += η*randn(seed, SVector{3})
            end
        else correlated
            a = arma_noise(N+1)
            for i in 1:length(tr)
                tr[i] = tr[i] .+ η .* a[i]
            end
        end
    end
    return tr
end
function rule_towel_noisy(x, p, n)
    @inbounds x1, x2, x3 = x[1], x[2], x[3]
    η, seed = p
    SVector(
        3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
        0.1*( (x2+0.35)*(1-2*x3) - 1 + η*randn(seed) )*(1 - 1.9*x1),
        3.78*x3*(1-x3) + 0.2*x2
    )
end
function towel_dynamic_noise(; N = default_N, η = 0.1, kwargs...)
    p = (η, MersenneTwister(919009))
    ds = DeterministicIteratedMap(rule_towel_noisy, [0.085, -0.121, 0.075], p)
    tr, = trajectory(ds, N; Ttr = 100)
    tr = standardize(tr)
end
function arma_noise(N)
    η = randn(N)
    s = ones(N)
    for n in 4:N
        s[n] = 1.625s[n-1] - 0.284s[n-2] - 0.355s[n-3] + η[n] - 0.96η[n-1]
    end
    s ./= std(s)
    return s
end

# Lorenz96
struct Lorenz96{N} end # Structure for size type
@inbounds function (obj::Lorenz96{N})(dx, x, p, t) where {N}
    dx[1] = (x[2] - x[N - 1]) * x[N] - x[1] + p[1]
    dx[2] = (x[3] - x[N]) * x[1] - x[2] + p[1]
    dx[N] = (x[1] - x[N - 2]) * x[N - 1] - x[N] + p[1]
    for n in 3:(N - 1)
        dx[n] = (x[n + 1] - x[n - 2]) * x[n - 1] - x[n] + p[1]
    end
    return nothing
end
function lorenz96_chaotic(; N = default_N, Δt = 0.2, D = 4, F = 24.0, kwargs...)
    ds = CoupledODEs(Lorenz96{D}(), range(0; length = D, step = 0.1), [F]; diffeq)
    tr, = trajectory(ds, N*Δt; Δt, Ttr = 1000.0)
    standardize(tr)
end
function lorenz96_rounded(; digits = 2, kwargs...)
    X = lorenz96_chaotic(; kwargs...)
    Z = [round.(x; digits) for x in X]
    return StateSpaceSet(Z)
end
function lorenz96_nonstationary_f!(dx, x, p, t)
    dF = p[1]
    N = length(x)-1
    F = x[end]
    @inbounds dx[1] = (x[2] - x[N - 1]) * x[N] - x[1] + F
    @inbounds dx[2] = (x[3] - x[N]) * x[1] - x[2] + F
    @inbounds dx[N] = (x[1] - x[N - 2]) * x[N - 1] - x[N] + F
    for n in 3:(N - 1)
      @inbounds dx[n] = (x[n + 1] - x[n - 2]) * x[n - 1] - x[n] + F
    end
    dx[end] = dF
    return nothing
end
function lorenz96_nonstationary(; N = default_N, Δt = 0.2, D = 6, F2=24.0, F1=1.0, kwargs...)
    u0 = range(0.1; length = D, step = 0.1)
    dF = (F2-F1)/(N*Δt)
    u0 = [u0..., F1]
    ds = CoupledODEs(lorenz96_nonstationary_f!, u0, [dF]; diffeq)
    tr, = trajectory(ds, N*Δt; Δt, Ttr = 0)
    standardize(tr)
end
function lorenz96_nonstationary_2(; N = default_N, Δt = 0.2, D = 6, F2=40.0, F1=10.0, kwargs...)
    tr1 = lorenz96_chaotic(; D, F = F1, N = N÷2, Δt)
    tr2 = lorenz96_chaotic(; D, F = F2, N = N÷2, Δt)
    return standardize(append!(tr1, tr2))
end

# Logistic maps
logistic(u, r) = r*u*(1-u)
function logistics!(un, u, p, t)
    L = length(u)
    k, r, η = p # η is dynamic noise
    for i in 2:L-1
        un[i] = logistic(u[i] + k*(u[i-1] - 2u[i] + u[i+1]) + η*randn(), r)
    end
    # circling around
    un[1] = logistic(u[1] + k*(u[L] - 2u[1] + u[2]) + η*randn(), r)
    un[L] = logistic(u[L] + k*(u[L-1] - 2u[L] + u[1]) + η*randn(), r)
    return
end
function coupled_logistics(; N = default_N, D = 8, k = 0.1, r = 4.0, η = 0, kwargs...)
    u0 = range(0.1, 0.2; length = D) .+ 0.0025671
    ds = DeterministicIteratedMap(logistics!, u0, [k, r, η])
    tr, = trajectory(ds, N; Ttr = 10000)
    return standardize(tr)
end

# Experimental
using DelimitedFiles
function experimental_data(; N = nothing, name, kwargs...)
    file = datadir("experimental", name*".txt")
    if name == "electrochemical1"
        τ = (0, 26, 13, 5, 20)
        x = vec(readdlm(file))
        A = genembed(x, τ)
    elseif name == "electrochemical2"
        τ = (0, 25, 16, 148, 138, 87, 60, 105)
        x = vec(readdlm(file))
        A = genembed(x, τ)
    elseif name == "shinriki1"
        τ = (0, 19, 38, 57)
        # js = ones(length(τ)) # perhaps use multivariate here
        x = readdlm(file)[:, 1]
        A = genembed(x, τ)
    elseif name == "shinriki2"
        τ = (0, 20, 40, 60)
        # js = ones(length(τ)) # perhaps use multivariate here
        x = readdlm(file)[:, 1]
        A = genembed(x, τ)
    elseif name == "nifty50"
        τ = [j*43 for j in 0:5]
        x = vec(readdlm(file))
        A = genembed(x, τ)
    elseif name == "vostok"
        τ = [j*50 for j in 0:7]
        x = vec(readdlm(file))
        A = genembed(x, τ)
    elseif name == "roessler_embed"
        τ = [0, 6, 3, 14]      # from PECUZAL
        τ = [j*7 for j in 0:6] # from TRADITIONAL
        x = readdlm(file)[:, 1]
        A = genembed(x, τ)
    elseif name == "roessler_all"
        M = readdlm(file)
        A = StateSpaceSet(M) # just take all 28 oscillators
    elseif name == "doublependulum"
        M = readdlm(file)
        τx = [0, 51, 25, 39, 12]
        τy = [0, 60, 30, 46, 14]
        jx = fill(1, length(τx))
        jy = fill(2, length(τy))
        τd = [0, 61, 31, 47, 31, 0]
        jd = [2, 2, 2, 2, 1, 1]
        A = StateSpaceSet(M)
        genembed(A, τx, jx)
    else
        error("Experimental dataset $(name) is uknown")
    end
    A = standardize(A)
    if isnothing(N)
        return A
    else
        return StateSpaceSet(A[1:min(length(A), N)])
    end
end

function electronic_roessler_mean_field(;R, Y = 1, d = 4, kwargs...)
    tr = mean(readdlm(datadir("exp_pro","R$(R)_ST_10_$Y.dat")); dims = 2)
    x = standardize(vec(tr))
    τ = estimate_delay(x, "mi_min")
    return embed(x, d, τ)
end

function electronic_roessler(; R=1, Y = 1, dms = 1:28, kwargs...)
    tr = readdlm(datadir("exp_pro","R$(R)_ST_10_$Y.dat"))[:,dms]
    standardize(StateSpaceSet(tr))
end

function embed_system(; system, d = 4, τ = nothing, kwargs...)
    data_producing_function = getfield(Data, system)
    x = data_producing_function(; kwargs...)[:, 1]
    x = standardize(x)
    if isnothing(τ)
        τ = estimate_delay(x, "mi_min")
    end
    return embed(x, d, τ)
end

end
