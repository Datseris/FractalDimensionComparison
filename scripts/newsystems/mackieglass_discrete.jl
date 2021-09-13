using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))
using3D()


function mackieglass(xnew, x, p, n)
    τ = length(x) # total delay, in units of δt
    β, γ, n, δt = p
    # x[1] is latest x
    xnew[1] = δt*(β*x[L]/(1 + x[L]^n) - γ*x[1]) + x[1]
    xnew[2:τ] .= view(x, 1:τ-1)
end

L = 300
p0 = [2, 1.0, 9.65, 0.01]
u0 = 0.1rand(L)

ds = DiscreteDynamicalSystem(mackieglass, u0, p0)
tr = trajectory(ds, 200000; Ttr = 10000, dt = 10)
x, y = tr[:, 1], tr[:, L]
figure()
plot(x,y, marker = "o", ls = "None", alpha = 0.5)
#
# t = estimate_delay(x, "mi_min")
# E = delay_afnn(x, t, 1:10)
# plot(1:10, E)
