using DynamicalSystems, TimeseriesSurrogates

D = 7
lo = Systems.lorenz96(D, range(0; length = D, step = 0.1); F = 8.0)
X = regularize(trajectory(lo, 1000, dt = 0.1, Ttr = 100.0))

e = 10.0 .^ range(-4, 1, length = 22)
CX = correlationsum(X, e; w = 5)

figure()
le = log10.(e)
plot(le, log10.(CX))

i = findfirst(z -> z > 0, CX)

@show linear_region(le[i:end], log10.(CX)[i:end])

sg = surrogenerator(X, ShuffleDimensions())
for i in 1:10
    Z = sg()
    CZ = correlationsum(Z, e)
    plot(le, log10.(CZ), alpha = 0.5, lw = 0.5, color = "k", ls = ":")
    @show linear_region(le, log10.(CZ))[2]
end
xlabel("log(e)")
ylabel("log(C)")
