# %% Estimation of least amount of points to get dimension
using DrWatson
@quickactivate :FractalDimensionComparison
include(srcdir("style.jl"))
using SpecialFunctions

# Eq. (3) from Tsonis, A. A., Elsner, J. B., & Georgakakos, K. P. (1993).
# Journal of the Atmospheric Sciences, 50(15), 2549–2555.
Nmin(n, k, A = 0.05) = sqrt(2)*sqrt(gamma(n/2 + 1))/(A*n*log(k))^((n+2)/2) *
    ( (2(k-1)*gamma(n+2)) / (gamma(0.5)^2 * gamma(n+3/2)) )^(n/2) * (n+2)/2

fig = figure()

ns = 2:0.25:8

for A in (0.05, 0.1)
    plot(ns, Nmin.(ns, 4, A); label = "A=$A")
end
legend()
