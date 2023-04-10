using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

ds = Systems.lorenz()
integ = integrator(ds)

function trajectory(ds::CDS{IIP, S, D}, integ, T, u, dt, Ttr) where {IIP, S, D}
    t0 = ds.t0
    tvec = (t0+Ttr):dt:(T+t0+Ttr)
    D = isnothing(a) ? D : length(a)
    sol = Vector{SVector{D, eltype(S)}}(undef, length(tvec))
    step!(integ, Ttr)
    for (i, t) in enumerate(tvec)
        while t > integ.t
            step!(integ)
        end
        sol[i] = SVector{D, eltype(S)}(obtain_access(integ(t), a))
    end
    return Dataset(sol)
end
