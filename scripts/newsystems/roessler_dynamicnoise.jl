using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
using StochasticDiffEq

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

# Initial conditions
u₀ = [0.1, -0.2, 0.1]
# short tspan just for plot demonstrating the noisiness for different η
tspan = (0.,100.)

for η in [0.1, 2]

    p = [0.2, 0.2, 5.7, η]
    prob_sde_roessler = SDEProblem(roessler_dynoise, σ_roessler_dynoise, u₀, tspan, p)

    sol = solve(prob_sde_roessler, LambaEM())

    fig = figure()
    u = hcat(sol.u...)'
    plot3D(u[:,1], u[:,2], u[:,3])
    name = savename("roessler_dynoise", @dict(η), "png")
    wsave(plotsdir(name),fig)
end
