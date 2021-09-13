### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ ad4957b2-5292-11eb-1f64-c96793eb50a3
using DrWatson

# ╔═╡ b993ce0a-52a4-11eb-325d-9762114e31e8
using BenchmarkTools, DynamicalSystems, Neighborhood, Entropies, DataFrames, CSV

# ╔═╡ bbc53e10-53f9-11eb-1166-17295832370f
using Plots

# ╔═╡ 3660b24a-52a5-11eb-06a7-5b4bf9b5f359
md"# Setup"

# ╔═╡ 02804130-5295-11eb-2076-9dd3d942b26c
@quickactivate "FractalDimension"

# ╔═╡ f08543ee-52a4-11eb-1928-89a52d10ffdf
md"Import all needed functions for benchmarks from packages and by direct definition."

# ╔═╡ 00c03bce-52a5-11eb-1fa0-19ffb17b23dd
function _correlationsum_new(tree, X, ε, w, qf)
    theiler = Neighborhood.Theiler(w)
    searchtype = Neighborhood.WithinRange(maximum(ε))
    N = length(X)
    Cq = similar(ε, Float64)
    @inbounds for i in 1:N
        idxs = Neighborhood.isearch(tree, X[i], searchtype, theiler(i))
        isempty(idxs) && continue
        dist = [Neighborhood.evaluate(tree.metric, X[i], X[j]) for j in idxs]
        sort!(dist)
        for k in 1:length(ε)-1
            j = searchsortedfirst(dist, ε[k])
            Cq[k] += qf(j-1)
        end
        Cq[end] += qf(length(dist))
    end
    return Cq ./ (N*qf(N-w-1))
end

# ╔═╡ fa4e2bfc-52a4-11eb-088f-9bd8377cc691
function correlationsum_new(X, ε; q = 2, metric = Euclidean(), w = 0)
    ntype = WithinRange(maximum(ε))
    tree = KDTree(X)
    qf = q == 2 ? identity : x -> float(x)^(q-1)
    _correlationsum_new(tree, X, ε, w, qf)
end

# ╔═╡ 07c6fc50-52a5-11eb-37d6-e546a66a5c28
# Define a function to access data sets.
"""
	sample_data(data, N_sample)
Sample `N_sample` unique random points out of `data`.
"""
function sample_data(data, N_sample)
    N_data = length(data)
    possible_indices = 1:N_data
    indices = unique(rand(possible_indices, 2N_sample))[1:N_sample]
    data[indices] |> Dataset
end

# ╔═╡ 3f6eb602-52a5-11eb-1735-f9d1f28fbcd8
md"# Single Benchmark Functions
## Classic Correlation"

# ╔═╡ 11e7cb6a-52a5-11eb-16de-cd668a8473ff
function benchmark_correlation(data)
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = minimum(maxima($data) - minima($data))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        [correlationsum($data, ε) for ε in εs]
    end
    minimum(bm.times)
end

# ╔═╡ 5a1a167a-52a5-11eb-11a9-958f6cefa7b3
md"## New Correlation"

# ╔═╡ 162b224e-52a5-11eb-0351-955a17fd1524
function benchmark_new_correlation(data)
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = minimum(maxima($data) - minima($data))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        correlationsum_new($data, εs)
    end
    minimum(bm.times)
end

# ╔═╡ 68a515f2-52a5-11eb-29ec-310b27c74895
md"## Boxed Correlation"

# ╔═╡ 1a849636-52a5-11eb-11fe-61afc4a4b453
function benchmark_boxed_correlation(data::Dataset{D,T}, prismdim = D) where {D,T}
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = estimate_r0_buenoorovio($data, $prismdim)
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        boxed_correlationsum($data, εs, ε1, M = $prismdim)
    end
    minimum(bm.times)
end

# ╔═╡ 7da26766-52a5-11eb-2318-593d397abb4c
md"## Moltenos Generalized Entropy"

# ╔═╡ 1f5dbdb8-52a5-11eb-3814-33ca2ad86123
function benchmark_molteno(data)
    bm = @benchmark molteno_boxing($data)
    minimum(bm.times)
end

# ╔═╡ 8a690dba-52a5-11eb-0703-932f4a1ab648
md"## Standard Generalized Entropy"

# ╔═╡ 2384e768-52a5-11eb-38cb-15b49cb5fdcb
function benchmark_genentropy(data)
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = minimum(maxima($data) - minima($data))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
		[Entropies.genentropy($data, ε) for ε in εs]
	end
	minimum(bm.times)
end

# ╔═╡ 9663d82a-52a5-11eb-3fe8-3f753118ab60
md"# Actual Experiments
## Benchmark of Size
The values are spread logarithmically"

# ╔═╡ 289e88c6-52a5-11eb-113a-b78dd25b0154
"""
	benchmark_N(
		lower_number,
		upper_number,
		number_of_numbers;
		base = 10,
		feedback = false
	)
Benchmarks `number_of_numbers` values in the interval `base^lower_number`, `base^higher_number`. These values are distributed logarithmically. If `feedback` is `true` there is a text feedback.
"""
function benchmark_N(
	lower_number,
	upper_number,
	number_of_numbers;
	base = 10,
	feedback = false
)
	result = DataFrame(N = Int[], Model = [], Method = [], Time = [])
	numbers = floor.(Int, base .^ range(lower_number, upper_number, length = number_of_numbers))

	# Define Systems.
	henon_sys = Systems.henon(rand(2))
	lorenz_sys = Systems.lorenz(randn(3))
	N_traj = 5 * floor(base ^ upper_number)

	# Integrate trajectories.
	henon = trajectory(henon_sys, N_traj, Ttr = 100)
	lorenz = trajectory(lorenz_sys, N_traj * 0.1, dt = 0.1, Ttr = 10)

	trajectories = [(:henon, henon), (:lorenz, lorenz)]

	for (name, traj) in trajectories
		feedback && println("Now using the $name system.")
		for N in numbers
			feedback && println("Currently calculating with $N random data points.")
			data = sample_data(traj, N)
			feedback && println("Benchmarking the classic correlation…")
			push!(
				result,
				(N, name, "Correlation",benchmark_correlation(data))
			)
			feedback && println("Benchmarking the new correlation…")
			push!(
				result,
				(N, name, "New Correlation", benchmark_new_correlation(data))
			)
			feedback && println("Benchmarking the boxed correlation…")
			push!(
				result,
				(N, name, "Boxed Correlation", benchmark_boxed_correlation(data))
			)
			feedback && println("Benchmarking the generalized entropy with factor 2 boxing…")
			push!(
				result,
				(N, name, "Molteno", benchmark_molteno(data))
			)
			feedback && println("Benchmarking the generalized entropy with standard boxing…")
			push!(
				result,
				(N, name, "Generalized Dimension", benchmark_genentropy(data))
			)
		end
	end
	wsave(
		datadir(
			"benchmark",
			savename(
				"benchmark_number",
				Dict(
					:lowerN => lower_number,
					:upperN => upper_number,
					:NofN => number_of_numbers,
					:base => base,
				),
				"csv",
			),
		),
		result,
	)
	result
end

# ╔═╡ 0ab12ebc-52a6-11eb-300d-1f864d2fcdf3
md"## Benchmark of dimension"

# ╔═╡ 2f2d2cec-52a5-11eb-0748-0dc568cf81aa
"""
	benchmark_dimension(
		lower_dimension::Int = 4,
		upper_dimension::Int = 16,
		distance_dimension::Int = 1;
		N_sample = 10000,
		feedback = false,
	)
Benchmarks different correlation algorithms for dimensions with distance `distance_dimension` between `lower_dimension` and `upper_dimension`. These algorithms are compared for different distance_dimensions. If `feedback` is `true` a status is printed after each benchmark step.
"""
function benchmark_dimension(
	lower_dimension::Int = 4,
	upper_dimension::Int = 16,
	distance_dimension::Int = 1;
	N_sample = 10000,
	feedback = false,
)
	result = DataFrame(Dimension = Int[], Method = [], Time = [])
	dimensions = lower_dimension:distance_dimension:upper_dimension
	for dimension in dimensions
		feedback && println("The current dimension is $dimension.")
		# The current trajectroy is longer than the data sample taken from the trajectory.
		traj = trajectory(
			Systems.lorenz96(dimension, rand(dimension), F = 8.0),
			N_sample,
			Ttr = 100,
			dt = 0.1,
		)
		# Only 10000 unique data points are sampled from the 100000 trajectory points
		data = sample_data(traj, N_sample)
		# Benchmark the classic correlation.
		feedback && println("Now benchmarking the classic correlation…")
		push!(
			result,
			(dimension, "Classic Correlation", benchmark_correlation(data)),
		)
		# Benchmark the new correlation that uses KDTrees.
		feedback && println("Now benchmarking the new correlation…")
		push!(
			result,
			(dimension, "New Correlation", benchmark_new_correlation(data)),
		)
		# Benchmark the Boxed Correlation without accessing the prism property.
		feedback && println("Now benchmarking the boxed correlation…")
		push!(
			result,
			(dimension, "Boxed Correlation", benchmark_boxed_correlation(data)),
		)
		feedback && println("Now benchmarking the prism correlation…")
		# Benchmark the Prism Correlation with Prism Dimension by Theiler.
		embedding_dim =
			0.75log2(N_sample) < dimension ?
			(0.5log2(N_sample) |> ceil |> Int) :
			dimension
		push!(
			result,
			(dimension, "Prism Correlation", benchmark_boxed_correlation(data, embedding_dim)),
		)
		feedback && println("Now benchmarking the generalized entropy by Molteno…")
		# Benchmark the generalized entropy method by Molteno.
		push!(
			result,
			(dimension, "Molteno", benchmark_molteno(data))
		)
		feedback && println("Now benchmarking the classic generalized entropy…")
		# Benchmark the classic generalized entropy.
		push!(
			result,
			(dimension, "Generalized Entropy", benchmark_genentropy(data))
		)
	end
	wsave(
		datadir(
			"benchmark",
			savename(
				"benchmark_dim",
				Dict(
					:N => N_sample,
					:lowerdim => lower_dimension,
					:upperdim => upper_dimension,
					:ddim => distance_dimension,
				),
				"csv",
			),
		),
		result,
	)
	result
end

# ╔═╡ fa38cdc0-52c0-11eb-08ba-e3b2e399dc2f
data = trajectory(Systems.henon(0.1rand(2)), 1000000, Ttr = 10)

# ╔═╡ fa33702a-52c0-11eb-250b-9d7a7a7506ad
benchmark_correlation(sample_data(data, 1000))

# ╔═╡ 2a048274-52c1-11eb-339c-6b10fb972d11
benchmark_new_correlation(sample_data(data, 1000))

# ╔═╡ 33866012-52c1-11eb-329c-eff5928eb12e
benchmark_boxed_correlation(sample_data(data, 1000))

# ╔═╡ 40b6a85a-52c1-11eb-080a-272e1642425d
benchmark_molteno(sample_data(data, 1000))

# ╔═╡ 4c2583a0-52c1-11eb-35ec-49b3c59c471f
benchmark_genentropy(sample_data(data, 1000))

# ╔═╡ 5c2b6d64-52c1-11eb-14e5-bd81d666bb81
benchmark_boxed_correlation(sample_data(data, 1000), 1)

# ╔═╡ 0124b0f2-53f9-11eb-3ebb-f182fd11a129
md"# Tests"

# ╔═╡ 64d806be-53f9-11eb-2ca9-1d2ece3404bd
sdata = sample_data(data, 10000)

# ╔═╡ e7bf3d5c-5414-11eb-0556-452471c2b65b
q = 1.99

# ╔═╡ bd059f48-5414-11eb-08c2-794f627f15ee


# ╔═╡ bba5b5d8-5415-11eb-35b1-c3e070d1e09d
md"# Results
## `boxed_correlationsum` vs `correlationsum`
It seems like they behave the same for a given ε. 

Say `henon` is a dataset with a trajectory of the henon map."

# ╔═╡ c6d939e8-5415-11eb-0c37-c73be604f29c
henon = trajectory(Systems.henon(0.1rand(2)), 10000, Ttr = 22)

# ╔═╡ 0d016116-5416-11eb-3ab5-a534b0d1453a
md"Then we can choose the εs by selecting the upper border `r0` using `estimate_r0_buenoorovio` and the lower border `ε_min` via `min_pairwise_distance` and spread them logarithmically on that interval."

# ╔═╡ 6f5e0968-5416-11eb-0fcb-572aaa902574
r0 = estimate_r0_buenoorovio(henon)

# ╔═╡ a0f41bf2-5416-11eb-2da6-cb4b485084c5
ε_min = min_pairwise_distance(henon)[2]

# ╔═╡ a2ca2444-5416-11eb-0a40-1f65148b816e
εs = 10 .^ range(log10(ε_min), log10(r0), length = 15)

# ╔═╡ 0798e090-53f9-11eb-377a-fdf6c2ade4d1
bcs = boxed_correlationsum(sdata, εs, q = q)

# ╔═╡ 88515afc-5415-11eb-273c-83e703c7392d
cs = correlationsum(sdata, εs, q = q)

# ╔═╡ 9fc5c036-5415-11eb-1d0e-5377b9d92943
bcs ≈ cs

# ╔═╡ be07006c-53f9-11eb-0141-4f9afae62d4b
plot(εs .|> log10, bcs .|> log10)

# ╔═╡ e89063f6-53fa-11eb-2671-5bacb746b401
linear_region(εs .|> log10, bcs .|> log10)[2]*1/(q-1)

# ╔═╡ efe5247c-5416-11eb-235e-6d7ce045b31e
md"Now the two correlation sums can be calculated."

# ╔═╡ fc2bd50a-5416-11eb-38f5-97335565a8eb
correlation_2 = correlationsum(henon, εs)

# ╔═╡ 067a4be0-5417-11eb-2774-915824fb7987
boxed_correlation_2 = boxed_correlationsum(henon, εs)

# ╔═╡ cbcac542-5416-11eb-3b99-791d68ce0745
md"Comparing these two directly works."

# ╔═╡ 1b3d11c0-5417-11eb-0856-b5aa0be46b0b
correlation_2 == boxed_correlation_2

# ╔═╡ 3842c2ec-5417-11eb-0030-83c97d87e157
md"The algorithm differs for `q≠2` so we need to compare this case as well."

# ╔═╡ 4d1e7eea-5417-11eb-1021-d50ee8615148
correlation_q = correlationsum(henon, εs, q = 3.3)

# ╔═╡ 61abb6b6-5417-11eb-32a7-f74a2549ef42
boxed_correlation_q = boxed_correlationsum(henon, εs, q = 3.3)

# ╔═╡ 70a0118a-5417-11eb-12cf-016cdbccb50d
correlation_q == boxed_correlation_q

# ╔═╡ 81510afc-5417-11eb-3111-6dfc809e8f6b
md"This does not come unexpected. The sums are only expected to be approximately equal, since the precision of the types plays a role."

# ╔═╡ c8c092ea-5417-11eb-2175-b3a1c62f1454
correlation_q ≈ boxed_correlation_q

# ╔═╡ 8f48f42e-5419-11eb-1865-c19b60d19b02
begin 
	import ChaosTools: boxed_correlationdim
	function boxed_correlationdim(data, εs, r0 = maximum(εs); q = 2, M = size(data, 2))
		@assert M ≤ size(data,2) "Prism dimension has to be lower or equal than " *
		"data dimension."
		if isinf(q)
			dd = boxed_correlationsum(data, εs, r0; q = 3, M = M)
			return 1 / (3 - 1) * linear_region(log.(εs), log.(maximum.(dd)), tol = 0.1)[2]
		else
			dd = boxed_correlationsum(data, εs, r0; q = q, M = M)
			return 1 / (q - 1) * linear_region(log.(εs), log.(dd), tol = 0.1)[2]
		end
	end
end

# ╔═╡ 420e26b2-5418-11eb-0477-5fd4284ff1d0
boxed_correlationdim(henon, q = 3)

# ╔═╡ e3627ac8-5417-11eb-1cf7-e13f24c1b62f
md"## q-order correlation literature
I was not able to find any paper using the explicit q-order correlationsum, yet its definition seems sensible and the plots I retrieve show linear scaling regions for different `q`."

# ╔═╡ d0e25b7c-5424-11eb-04dd-396f8e034b86
begin 
	ϵ_min = min_pairwise_distance(henon)[2]
	ϵ_max = maximum(maxima(henon)-minima(henon))
	ϵs = 10 .^ range(log10(ϵ_min), log10(ϵ_max), length = 20)
	plot()
	for q in 1.5:0.5:6
		corsum = correlationsum(henon, ϵs, q = q)
		(p1, p2), dim = linear_region(log10.(ϵs), log10.(corsum))
		dim = 1/(q-1)*dim
		plot!(
			log10.(ϵs), 
			log10.(corsum), 
			label = "q = $q, Δ = $(round(dim, digits = 3))",
		)
	end
	plot!(
		legend = :bottomright, 
		xlims = log10.((ϵs[1], ϵs[20])), 
		xguide = "\$ \\log(\\epsilon) \$",
		yguide = "\$ C_q \$",
	)
end
		

# ╔═╡ Cell order:
# ╟─3660b24a-52a5-11eb-06a7-5b4bf9b5f359
# ╠═ad4957b2-5292-11eb-1f64-c96793eb50a3
# ╠═02804130-5295-11eb-2076-9dd3d942b26c
# ╟─f08543ee-52a4-11eb-1928-89a52d10ffdf
# ╠═b993ce0a-52a4-11eb-325d-9762114e31e8
# ╠═fa4e2bfc-52a4-11eb-088f-9bd8377cc691
# ╠═00c03bce-52a5-11eb-1fa0-19ffb17b23dd
# ╠═07c6fc50-52a5-11eb-37d6-e546a66a5c28
# ╟─3f6eb602-52a5-11eb-1735-f9d1f28fbcd8
# ╠═11e7cb6a-52a5-11eb-16de-cd668a8473ff
# ╟─5a1a167a-52a5-11eb-11a9-958f6cefa7b3
# ╠═162b224e-52a5-11eb-0351-955a17fd1524
# ╟─68a515f2-52a5-11eb-29ec-310b27c74895
# ╠═1a849636-52a5-11eb-11fe-61afc4a4b453
# ╟─7da26766-52a5-11eb-2318-593d397abb4c
# ╠═1f5dbdb8-52a5-11eb-3814-33ca2ad86123
# ╟─8a690dba-52a5-11eb-0703-932f4a1ab648
# ╠═2384e768-52a5-11eb-38cb-15b49cb5fdcb
# ╟─9663d82a-52a5-11eb-3fe8-3f753118ab60
# ╠═289e88c6-52a5-11eb-113a-b78dd25b0154
# ╟─0ab12ebc-52a6-11eb-300d-1f864d2fcdf3
# ╠═2f2d2cec-52a5-11eb-0748-0dc568cf81aa
# ╠═fa38cdc0-52c0-11eb-08ba-e3b2e399dc2f
# ╠═fa33702a-52c0-11eb-250b-9d7a7a7506ad
# ╠═2a048274-52c1-11eb-339c-6b10fb972d11
# ╠═33866012-52c1-11eb-329c-eff5928eb12e
# ╠═40b6a85a-52c1-11eb-080a-272e1642425d
# ╠═4c2583a0-52c1-11eb-35ec-49b3c59c471f
# ╠═5c2b6d64-52c1-11eb-14e5-bd81d666bb81
# ╟─0124b0f2-53f9-11eb-3ebb-f182fd11a129
# ╠═64d806be-53f9-11eb-2ca9-1d2ece3404bd
# ╠═e7bf3d5c-5414-11eb-0556-452471c2b65b
# ╠═0798e090-53f9-11eb-377a-fdf6c2ade4d1
# ╠═88515afc-5415-11eb-273c-83e703c7392d
# ╠═9fc5c036-5415-11eb-1d0e-5377b9d92943
# ╠═bbc53e10-53f9-11eb-1166-17295832370f
# ╠═be07006c-53f9-11eb-0141-4f9afae62d4b
# ╠═e89063f6-53fa-11eb-2671-5bacb746b401
# ╠═bd059f48-5414-11eb-08c2-794f627f15ee
# ╟─bba5b5d8-5415-11eb-35b1-c3e070d1e09d
# ╠═c6d939e8-5415-11eb-0c37-c73be604f29c
# ╟─0d016116-5416-11eb-3ab5-a534b0d1453a
# ╠═6f5e0968-5416-11eb-0fcb-572aaa902574
# ╠═a0f41bf2-5416-11eb-2da6-cb4b485084c5
# ╠═a2ca2444-5416-11eb-0a40-1f65148b816e
# ╟─efe5247c-5416-11eb-235e-6d7ce045b31e
# ╠═fc2bd50a-5416-11eb-38f5-97335565a8eb
# ╠═067a4be0-5417-11eb-2774-915824fb7987
# ╟─cbcac542-5416-11eb-3b99-791d68ce0745
# ╠═1b3d11c0-5417-11eb-0856-b5aa0be46b0b
# ╟─3842c2ec-5417-11eb-0030-83c97d87e157
# ╠═4d1e7eea-5417-11eb-1021-d50ee8615148
# ╠═61abb6b6-5417-11eb-32a7-f74a2549ef42
# ╠═70a0118a-5417-11eb-12cf-016cdbccb50d
# ╟─81510afc-5417-11eb-3111-6dfc809e8f6b
# ╠═c8c092ea-5417-11eb-2175-b3a1c62f1454
# ╠═8f48f42e-5419-11eb-1865-c19b60d19b02
# ╠═420e26b2-5418-11eb-0477-5fd4284ff1d0
# ╟─e3627ac8-5417-11eb-1cf7-e13f24c1b62f
# ╟─d0e25b7c-5424-11eb-04dd-396f8e034b86
