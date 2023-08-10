# %% Estimation of least amount of points to get dimension
using DrWatson
@quickactivate :FractalDimensionComparison

using SpecialFunctions

# Eq. (3) from Tsonis, A. A., Elsner, J. B., & Georgakakos, K. P. (1993).
# Journal of the Atmospheric Sciences, 50(15), 2549–2555.
Nmin(n, k, A = 0.05) = sqrt(2)*sqrt(gamma(n/2 + 1))/(A*n*log(k))^((n+2)/2) *
    ( (2(k-1)*gamma(n+2)) / (gamma(0.5)^2 * gamma(n+3/2)) )^(n/2) * (n+2)/2

fig = Figure()
ax = Axis(fig[1,1])
ns = 2:0.25:8

for A in (0.05, 0.1)
    lines!(ns, Nmin.(ns, 4, A); label = "A=$A")
end
axislegend(ax)

# %% Diameters of KS, L32
εmax = MathConstants.e^2
using LinearAlgebra

X = Data.lorenz96_chaotic(D=32, F = 24.0)
mini, maxi = minmaxima(X)
dL32 = norm(maxi .- mini)

X = Data.ksiva()
mini, maxi = minmaxima(X)
dKS = norm(maxi .- mini)
