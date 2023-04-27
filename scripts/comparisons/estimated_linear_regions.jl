using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff
include(srcdir("style.jl"))

datasets = Vector(undef, 6)

X = Data.roessler_chaotic(; N = 1e4)
εs = MathConstants.e .^ (-8:0.2:5)

Hs = genentropy.(Ref(X), εs)

x = log.(εs)
y = - Hs
figure()
plot(x, y)

lrs, slopes = linear_regions(x, y; tol = 0.25)

plot(x, y)
for i in 1:length(lrs)-1
    plot(
        x[lrs[i]:lrs[i+1]], y[lrs[i]:lrs[i+1]];
        marker = "o", ms = 5
    )
end