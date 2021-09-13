#=
This file attempts to implement the "imporoved dimension estimation" by Judd:

Judd, K. (1992). An improved estimator of dimension and some comments on providing
confidence intervals. Physica D: Nonlinear Phenomena, 56(2–3), 216–228. https://doi.org/10.1016/0167-2789(92)90025-I

Judd, K. (1994). Estimating dimension from small samples. Physica D: Nonlinear Phenomena,
71(4), 421–429. https://doi.org/10.1016/0167-2789(94)90008-6

So far this attempt has failed.
=#
using Distances, BigCombinatorics


"""
    interpoint_distances(X, N, norm = Euclidean()) -> ds

returns N interpoint distances of randomly chosen pairs of points from time series X.
N must be smaller than the free memory. N ≈ 5000 is actually already sufficient for the histograms. N ≈ 45000 sort of the upper limit of what you should use, this makes for a smoother hist later
"""
function interpoint_distances(X, N, norm = Euclidean())
    ds = eltype(X)[]

    @assert N < Int(Sys.free_memory()) ÷ 8

    # choosing N random points from X.
    sizehint!(ds, Int(N*(N-1)/2))
    while length(ds) < Int(N*(N-1)/2)
        i,j = rand(1:length(X),2)
        dist = evaluate(norm, X[i], X[j])
        (i != j && !any(ds .== dist)) && push!(ds, dist)
    end
    return filter(x-> x!= 0.0,ds)
end

"""
    logspace_histogram(distances::Vector, λ = exp(-w0), ϵ₀ = λ^2 * maximum(x)) -> bins, counts
Calculate a logarithmically spaced histogram of the values in `x` (must be positive
definite, i.e. distances) and return the bins and count of values in the bins.

The `i`-th bin is `[ϵ₀*λ^i, ϵ₀*λ^(i+1))`. Elements of `x` greater than `ϵ₀` are disregarded.
"""
function logspace_histogram(distances::AbstractVector, ε₀, λ = exp(-w0(x)))
    @assert λ < 1
    εs = [ε₀]
    λmin = minimum(distances)*λ
    εᵢ = ε₀
    while εᵢ > λmin
        εᵢ = εs[end]*λ
        push!(εs, εᵢ)
    end

    reverse!(εs)
    L = length(εs)
    bs = zeros(Int, L)
    # Fill bins with values
    for v in distances
        i = searchsortedfirst(εs, v)
        i == L + 1 && continue
        bs[i] += 1
    end
    return reverse!(εs), reverse!(bs)
end
# TODO: test the algorithm with including the bin [∞, ε0)
w0(x) = log(maximum(x)/minimum(x))/√length(x)


"""
    bin_judd(x::AbstractVector) -> εs, bs, ε₀

returns a logspace histogram with optimal bin width.
Estimation of the latter could probably be done more efficient

"""
function bin_judd(x::AbstractVector)

    opt_λ =  optimal_λ(x) # This is probably not very efficient
    εs, bs = logspace_histogram(x,maximum(x),opt_λ)

    ε₀ = εs[findmax(bs)[2]+1]
    new_εs = reverse(εs[findmax(bs)[2]+1:end])
    new_bs = reverse(bs[findmax(bs)[2]+1:end])
    push!(new_bs, sum(bs[1:findmax(bs)[2]]))
    push!(new_εs, εs[findmax(bs)[2]])
    reverse!(new_bs)
    reverse!(new_εs)
    return new_εs, new_bs, ε₀
end

function P(εᵢ, d,a,ε₀)
    polynomial = length(a) == 0 ? 1 : 0
    for t in 1:length(a)
        polynomial += a[t] * (εᵢ/ε₀)^(t-1)
    end
    return (εᵢ/ε₀)^d * polynomial
end

""" log_l(εs, bs, P, d, a,ε₀)

returns the log-likelihood function according to Judd
"""
log_l(εs, bs, P, d, a,ε₀) = -sum(bs[i]*log( P(εs[i],d,a,ε₀) - P(εs[i+1],d,a,ε₀) ) for i in 1:length(εs)-1 )

""" optimal_λ(dists)

searches for the optimal bin width w = log(1/λ) of the judd histogram.
The search is currently not how Judd proposed, because I didn't get anywhere with that. Right now, this is a brute force search, this should definitely be improved if this is actuallly going to be used.

"""
function optimal_λ(dists)
    n = length(dists)
    w₀ = w0(dists)
    ws = collect(w₀:0.001:.5)
    f_ws = []

    for w in ws
        εs, bs = logspace_histogram(dists, maximum(dists), exp(-w))
        m = findlast(x->x!=0,bs)
        f = log(w) + log(Multinomial(bs[1:m])) + log(Binomial(n+m+1,m))
        push!(f_ws,f)
    end

    w = ws[findmin(f_ws)[2]]
    λ = exp(-w)
    return λ
end

problargerzero(εs, bs, d₀,a₀,ε₀) = all(P(εs[i], d₀, a₀,ε₀)-P(εs[i+1],d₀,a₀,ε₀) > 0 for i in 1:length(εs)-1)

""" nested_optimizer(εs, bs; ε₀, d₀, a₀, stp = 1e-5, η = 1e-3, tol=1e-4, decay=1e-4, maxiter = 10000)

    uses numerical (center-diff) gradient descent with adaptive (decaying) learning rate to find
    a minimum of the log-liklelihood function.
    Seems to produce reasonable results if `deg(a) < 2`, but still sensitive to initial
    conditions.
    The key was to use the nested version that Judd suggested in a footnote of the 1994 paper:
    nest the optimizations so that first `d` is optimized with fixed `a`, then the coefficients
    of `a` are optimized with fixed `d`. Repeat until total step size is smaller than `tol`.

    There are a lot of parameters to this optimization problem that themselves can probably be
    optimized and/ or made more adaptive so that they fit the specific problem.

"""
function nested_optimizer(εs, bs; ε₀, d₀, a₀, stp = 1e-5, η = 1e-3, γ = 1e-3, tol=1e-4, decay=1e-4, maxiter = 10000)

    iter = 0

    deg_a = length(a₀) # Technically, this is not the *degree*, but deg(a)+1, but it's more convenient to write it like this

    cur_x = [d₀, a₀...]
    cur_d = cur_x[1]
    cur_a = cur_x[2:end]

    step_size = 10

    while iter < maxiter && step_size > tol

        prev_v = cur_x

        # First we will optimize d with fixed polynomial coefficients
        cur_d = cur_x[1]
        v = 0
        iter_d = 0
        d_step = 1
        while iter_d < maxiter÷10 && d_step > tol/10
            prev_d = cur_d
            g = (log_l(εs, bs, P, cur_d+stp, a₀,ε₀)-log_l(εs, bs, P, cur_d-stp, a₀,ε₀))/(2stp)
            v = η * g + γ * v
            cur_d -= v
            d_step = abs(cur_d-prev_d)
            iter_d += 1
        end

        # now we will optimize the coefficients of a with the previous best d
        grad = zeros(deg_a)
        cur_a = cur_x[2:end]
        v = zeros(deg_a)
        iter_a = 0
        a_step = 1
        while iter_a < maxiter÷10 && a_step > tol/10
            prev_a = cur_a
            for k in 1:deg_a
                a_stp = zeros(deg_a)
                a_stp[k] = stp
                grad[k] = (log_l(εs, bs, P, d₀, a₀ + a_stp,ε₀)-log_l(εs, bs, P, d₀, a₀-a_stp,ε₀))/(2stp)
            end
            v = γ * v + η * grad
            cur_a -= v
            a_step = evaluate(Euclidean(), cur_a, prev_a)
            iter_a += 1
        end

        if problargerzero(εs, bs, cur_d,cur_a,ε₀)
            cur_x = [cur_d, cur_a...]
        else
            break
        end

        # check total step size of parameter vector
        step_size = evaluate(Euclidean(), cur_x, prev_v)

        # update learning rate
        η /= (1+decay*iter)

        iter += 1
    end
    cur_d, cur_a
end
