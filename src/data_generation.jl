module Data

using DrWatson
using Random, Statistics
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra
using DelimitedFiles, DrWatson
using StochasticDiffEq

# default constants used if no other is given
const default_N = 10000
const diffeq = (alg = Vern9(), reltol = 1e-12, abstol = 1e-12, maxiters = typemax(Int))
const diffeq_lowacc = (alg = Tsit5(), reltol = 1e-9, abstol = 1e-9, maxiters = typemax(Int))

function kaplanyorke_map(; N = default_N, λ = 0.2, kwargs...)
    f(u, λ, t) = SVector((2u[1]) % 1.0 + 1e-15rand(Random.GLOBAL_RNG), λ*u[2] + cospi(4u[1]))
    ds = DiscreteDynamicalSystem(f, SVector(0.15, 0.2), λ)
    tr = trajectory(ds, N; Ttr = 1000)
    standardize(tr)
end

function roessler_periodic(; N = default_N, Δt = 0.5, kwargs...)
    ds = Systems.roessler([0.1, -0.2, 0.1]; c = 3.0, a = 0.2, b = 0.2)
    tr = trajectory(ds, N*Δt; Δt, Ttr = 1000, diffeq...)
    standardize(tr)
end

function roessler_chaotic(; η=0, correlated=false, N = default_N, Δt = 0.5, kwargs...)
    ds = Systems.roessler([0.1, -0.2, 0.1]; c = 5.7, a = 0.2, b = 0.2)
    tr = trajectory(ds, N*Δt; Δt, Ttr = 1000, diffeq...)
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
    return Dataset(Z)
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
    X = Dataset(hcat(sol.u...)')
    return standardize(X)
end

function roessler_lorenz(; N = default_N, kwargs...)
    tr1 = roessler_periodic(N = N÷2)
    tr2 = lorenz96_chaotic(; D = 4, N = N÷2)[:, 1:3]
    tr2 = trajectory(ds, (N/2)*Δt; Δt, Ttr = 1000, diffeq...)
    append!(standardize(tr1), standardize(tr2))
end

function roessler_sphere(; N = default_N, Δt = 0.5, kwargs...)
    ds = Systems.roessler([0.1, -0.2, 0.1]; c = 3.0, a = 0.2, b = 0.2)
    tr1 = standardize(trajectory(ds, (N/2)*Δt; Δt, Ttr = 1000, diffeq...))
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
    ds = ContinuousDynamicalSystem(roessler_nonstationary_eom, u0, p0)
    tr = trajectory(ds, N*Δt; Δt, Ttr = 0, diffeq...)
end

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
    return Dataset([torus(u) for u in zip(θs, φs)])
end

function koch(; maxk = 7, kwargs...)
    flakepoints = SVector{2}.([[0.0; 0.0], [0.5; sqrt(3)/2], [1; 0.0], [0.0; 0.0]])
    function innerkoch(points, maxk, α = sqrt(3)/2)
        Q = @SMatrix [0 -1; 1 0]
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
    kochdata = Dataset(unique!(kochpoints))
end

function henonheiles_quasi(; N = default_N, Δt = 0.2, kwargs...)
    u0 = [0.0, 0.1, 0.5, 0.0]
    ds = Systems.henonheiles(u0)
    tr = trajectory(ds, N*Δt; Δt, Ttr = 100, diffeq...)
    standardize(tr)
end

function henonheiles_chaotic(; N = default_N, Δt = 0.5, kwargs...)
    u0 = [0.0, -0.25, 0.42081, 0.0]
    ds = Systems.henonheiles(u0)
    tr = trajectory(ds, N*Δt; Δt, Ttr = 100, diffeq...)
    standardize(tr)
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
    return Dataset(A)
end

function standardmap_chaotic(; N = default_N, k = 64.0, kwargs...)
    ds = Systems.standardmap([0.08152, 0.122717]; k)
    tr = trajectory(ds, N; Ttr = 100)
end

function henon_chaotic(; N = default_N, a = 1.4, b = 0.3, kwargs...)
    ds = Systems.henon([0.08152, 0.122717]; a, b)
    tr = trajectory(ds, N; Ttr = 100)
    standardize(tr)
end

function towel_chaotic(; N = default_N, η = 0, correlated = false, kwargs...)
    ds = Systems.towel([0.085, -0.121, 0.075])
    tr = trajectory(ds, N; Ttr = 100)
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

function arma_noise(N)
    η = randn(N)
    s = ones(N)
    for n in 4:N
        s[n] = 1.625s[n-1] - 0.284s[n-2] - 0.355s[n-3] + η[n] - 0.96η[n-1]
    end
    s ./= std(s)
    return s
end


function rule_towel_noisy(x, p, n)
    @inbounds x1, x2, x3 = x[1], x[2], x[3]
    η, seed = p
    SVector(
        3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
        0.1*( (x2+0.35)*(1-2*x3) - 1 + η*randn(seed) )*(1 - 1.9*x1),
        3.78*x3*(1-x3)+0.2*x2
    )
end
function towel_dynamic_noise(; N = default_N, η = 0.1, kwargs...)
    p = [η, MersenneTwister(919009)]
    ds = DiscreteDynamicalSystem(rule_towel_noisy, [0.085, -0.121, 0.075], p)
    tr = trajectory(ds, N; Ttr = 100)
    tr = standardize(tr)
end


function lorenz96_chaotic(; N = default_N, Δt = 0.2, D = 4, F = 24.0, kwargs...)
    ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F)
    tr = trajectory(ds, N*Δt; Δt, Ttr = 1000.0, diffeq...)
    standardize(tr)
end

function lorenz96_rounded(; digits = 2, kwargs...)
    X = lorenz96_chaotic(; kwargs...)
    Z = [round.(x; digits) for x in X]
    return Dataset(Z)
end

function lorenz96_nonstationary_f!(dx, x, p, t)
    dF = p[1]
    N = length(x)-1
    F = x[end]

    # 3 edge cases
    @inbounds dx[1] = (x[2] - x[N - 1]) * x[N] - x[1] + F
    @inbounds dx[2] = (x[3] - x[N]) * x[1] - x[2] + F
    @inbounds dx[N] = (x[1] - x[N - 2]) * x[N - 1] - x[N] + F
    # then the general case
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
    ds = ContinuousDynamicalSystem(lorenz96_nonstationary_f!, u0, [dF])
    tr = trajectory(ds, N*Δt; Δt, Ttr = 0,  diffeq...)
    # standardize(tr)
end

function lorenz96_nonstationary_2(; N = default_N, Δt = 0.2, D = 6, F2=40.0, F1=10.0, kwargs...)
    ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = F1)
    tr1 = trajectory(ds, (N/2)*Δt; Δt, Ttr = 1000.0, diffeq...)

    ds = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = F2)
    tr2 = trajectory(ds, (N/2)*Δt; Δt, Ttr = 1000.0, diffeq...)

    return standardize(append!(tr1, tr2))
end

logistic(u) = 4*u*(1-u)
function logistics!(un, u, p, t)
    L = length(u)
    k, r, η = p # η is dynamic noise
    for i in 2:L-1
        un[i] = logistic(u[i] + k*(u[i-1] - 2u[i] + u[i+1]) + η*randn())
    end
    # circling around
    un[1] = logistic(u[1] + k*(u[L] - 2u[1] + u[2]) + η*randn())
    un[L] = logistic(u[L] + k*(u[L-1] - 2u[L] + u[1]) + η*randn())
    return
end

function coupled_logistics(; N = default_N, D = 8, k = 0.1, r = 4.0, η = 0, kwargs...)
    u0 = range(0.1, 0.9; length = D) .+ 0.0001
    ds = DiscreteDynamicalSystem(logistics!, u0, [k, r, η])
    tr = trajectory(ds, N; Ttr = 1000)
end

function embed_system(; system, d = 4, kwargs...)
    data_producing_function = getfield(Data, system)
    x = data_producing_function(; kwargs...)[:, 1]
    x = standardize(x)
    τ = estimate_delay(x, "mi_min")
    return embed(x, d, τ)
end

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
        A = Dataset(M) # just take all 28 oscillators
    elseif name == "doublependulum"
        M = readdlm(file)
        τx = [0, 51, 25, 39, 12]
        τy = [0, 60, 30, 46, 14]
        jx = fill(1, length(τx))
        jy = fill(2, length(τy))
        τd = [0, 61, 31, 47, 31, 0]
        jd = [2, 2, 2, 2, 1, 1]
        A = Dataset(M)
        genembed(A, τx, jx)
    else
        error("Experimental dataset $(name) is uknown")
    end
    A = standardize(A)
    if isnothing(N)
        return A
    else
        return Dataset(A[1:min(length(A), N)])
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
    standardize(Dataset(tr))
end

end
