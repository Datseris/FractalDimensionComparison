using DrWatson
@quickactivate :FractalDimension # uses DynamicalSystems, PyPlot
include(srcdir("style.jl"))
using3D()

function crazy_towel(x, p, n)
    x1, x2, x3, x4, x5 = x
    SVector(
        3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
        0.1*( (x2+0.35)*(1-2*x3) - 1 )*(1 - 1.9*x1),
        3.78*x3*(1-x3) + 0.5*x2*x4,
        2.0*x4*(1-x4*x3) - 0.2x3,
        0.2*( x5^2*(1-2*x1) - 1 )*(1 - 0.9*x4),
    )
end

u0 = -0.01rand(5)
p0 = rand(3)
ds = DiscreteDynamicalSystem(crazy_towel, u0, p0)

tr = trajectory(ds, 2000; Ttr = 10)
x, y, z, w, e = columns(tr)
scatter3D(z, w, e)
