using HypothesisTests, FractalDimensions
using PredefinedDynamicalSystems
using ComplexityMeasures
using CairoMakie

# ds = Systems.lorenz()
ds = Systems.lorenz96(6; F = 24.0)
X, t = trajectory(ds, 1000.0; Dt = 0.1, Ttr = 10.0)

p = 0.99
estimator = :mm
TestType = ExactOneSampleKSTest
sigmas, xis, pvalues = extremevaltheory_gpdfit_pvalues(A, p; TestType, estimator)

fig, ax = hist(1 ./ sigmas; axis = (title = "1 / sigmas",), bins = 50)
hist(fig[2,1], xis; axis = (title = "xis",), bins = 50)
hist(fig[1,2], pvalues; axis = (title = "pvalues",), bins = 50)
text(fig[2,2], "system: Lorenz96 (6 dim.)\np = $(p)\nestimator = $(estimator)\nTestType = $(TestType)"; align = (:center, :center))
fig

# %% Plot exemplary distributions
using FractalDimensions.ProgressMeter

ds = Systems.lorenz96(6; F = 24.0)
X, t = trajectory(ds, 1000.0; Dt = 0.1, Ttr = 10.0)

p = 0.95
estimator = :mm
TestType = ExactOneSampleKSTest

TestType = ApproximateOneSampleKSTest
N = length(X)
progress = ProgressMeter.Progress(
    N; desc = "Extreme value theory p-values: ", enabled = true
)
logdists = [zeros(eltype(X), N) for _ in 1:Threads.nthreads()]
pvalues = zeros(N)
nrmses = zeros(N)
sigmas = zeros(N)
xis = zeros(N)
Es = [Float64[] for _ in 1:N]

function gpd_nrmse(E, gpd; nbins = 20)

    bins = range(0, nextfloat(maximum(E), 4); length = nbins)
    binning = FixedRectangularBinning(bins)
    allprobs = allprobabilities(ValueHistogram(binning), E)

    # Calculate NRMSE. obtain mean value of each bin
    width = step(bins)
    midpoints = bins .+ width/2

    rmse = zero(eltype(E))
    mmse = zero(eltype(E))
    meandens = mean(allprobs)/width
    for (j, prob) in enumerate(allprobs)
        dens = prob / width # density value
        dens_gpd = pdf(gpd, midpoints[j])
        rmse += (dens - dens_gpd)^2
        mmse += (dens - meandens)^2
    end
    nrmse = sqrt(rmse/mmse)
    return nrmse
end

Threads.@threads for j in eachindex(X)
    logdist = logdists[Threads.threadid()]
    @inbounds map!(x -> -log(FractalDimensions.euclidean(x, X[j])), logdist, vec(X))
    σ, ξ, E = FractalDimensions.extremevaltheory_local_gpd_fit(logdist, p, estimator)
    sigmas[j] = σ
    xis[j] = ξ
    Es[j] = E
    # Note that exceedances are defined with 0 as their minimum
    gpd = FractalDimensions.GeneralizedPareto(0, σ, ξ)
    test = TestType(E, gpd)
    pvalues[j] = pvalue(test)
    nrmses[j] = gpd_nrmse(E, gpd)
    ProgressMeter.next!(progress)
end

# %% Analyze/plot/whatever
using Distributions

fig = Figure()

for i in 1:3
ax = Axis(fig[i, 1])
j = rand(eachindex(Es))
E = Es[j]
s, x, pval = sigmas[j], xis[j], pvalues[j]
gpd = GeneralizedPareto(0, s, x)
nbins = 25

hist!(ax, E; bins = nbins, normalization = :pdf)
xrange = range(0, maximum(E); length = 100)
nrmse = nrmses[j]
t = "pval = $(pval)\nnrmse = $(nrmse)"
lines!(ax, xrange, pdf.(gpd, xrange); label = t, color = :black)
axislegend(ax)
end

axh = Axis(fig[:, 2])
hist!(axh, pvalues; label = "pvalues", normalization = :pdf)
hist!(axh, nrmses; label = "nrmses", normalization = :pdf)
axislegend(axh)

Label(fig[0, :], "p = $(p)"; tellwidth = false)
display(fig)



# %%
using Distributions, HypothesisTests, CairoMakie

sigma = 1 / 2.0
xi = -0.1

gpd = GeneralizedPareto(0.0, sigma, xi)

X = rand(gpd, 10000)

TestType = ExactOneSampleKSTest
test = TestType(X, gpd)

p = pvalue(test)

fig, ax = hist(X; bins = 50, normalization = :pdf, label = "pvalue = $(round(p; digits=3))")
xrange = range(0, maximum(X); length = 100)
lines!(xrange, pdf.(gpd, xrange); color = :black, label = "analytic")
axislegend(ax)
display(fig)

# test Cramer direct implementation
n = length(X)
xs = sort(X)

T = 1/(12n) + sum(i -> (cdf(gpd, xs[i]) - (2i - 1)/2n)^2, 1:n)


# An approximation of the T statistic is normally distributed
# with std of approximately sqrt(1/45)
# From there the z statistic
zstat = T/sqrt(1/45)

newp = 2*(1 - cdf(Normal(0,1), zstat))
# # it is easy to approxamte the T statistic
# # omega2 is Chi Squared distributed
# omegadist = Chisq(n-1)
# omegaquant = quantile(omegadist, 0.95)
@show newp

# %% Test with normal

using Distributions, HypothesisTests

sigma = 1 / 2.0
m = 4.6

d = Normal(m, sigma)

X = rand(d, 1000)

TestType = OneSampleADTest
test = TestType(X, d)

p = pvalue(test)

fig, ax = hist(X; bins = 50, normalization = :pdf, label = "pvalue = $(round(p; digits=3))")
xrange = range(0, maximum(X); length = 100)
lines!(xrange, pdf.(d, xrange); color = :black, label = "analytic")
axislegend(ax)
fig


# %% Test that p-value is uniformly distributed if hypothesis is true

using HypothesisTests, Distributions

sigma = 1 / 2.0
xi = -0.1
n = 1000 # sample size

gpd = GeneralizedPareto(0.0, sigma, xi)

pvals = []

for i in 1:10000

    X = rand(gpd, 10000)

    TestType = ExactOneSampleKSTest
    test = TestType(X, gpd)

    p = pvalue(test)
    push!(pvals, p)
end

using CairoMakie

hist(pvals)