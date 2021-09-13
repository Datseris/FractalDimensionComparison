using DrWatson
@quickactivate "FractalDimension"

# Import all needed functions for benchmarks from packages and by direct definition.
using BenchmarkTools, DynamicalSystems, Neighborhood, Entropies, DataFrames, CSV

function correlationsum_new(X, ε; q = 2, metric = Euclidean(), w = 0)
    ntype = WithinRange(maximum(ε))
    tree = KDTree(X)
    qf = q == 2 ? identity : x -> float(x)^(q-1)
    _correlationsum_new(tree, X, ε, w, qf)
end

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

# Define Benchmark functions.
# All functions can output probabilities so the time is benchmarked for probabilities.
function benchmark_correlation(data)
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = minimum(maxima($data) - minima($data))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        [correlationsum($data, ε) for ε in εs]
    end
    minimum(bm.times)
end

function benchmark_new_correlation(data)
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = minimum(maxima($data) - minima($data))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        correlationsum_new($data, εs)
    end
    minimum(bm.times)
end

function benchmark_boxed_correlation(data::Dataset{D,T}, prismdim = D) where {D,T}
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = estimate_r0_buenoorovio($data, $prismdim)
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
        boxed_correlationsum($data, εs, ε1, M = $prismdim)
    end
    minimum(bm.times)
end

function benchmark_molteno(data)
    bm = @benchmark molteno_boxing($data)
    minimum(bm.times)
end

function benchmark_genentropy(data)
    bm = @benchmark begin
        ε0 = min_pairwise_distance($data)[2]
        ε1 = minimum(maxima($data) - minima($data))
        εs = 10 .^ range(log10(ε0), log10(ε1), length = 12)
		[Entropies.genentropy($data, ε) for ε in εs]
	end
	minimum(bm.times)
end

# Define size benchmark with logarithmically spread values.
# `base = 1` allows non-logarithmic values
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

# Define Dimension Benchmark.
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

# df_dim = benchmark_dimension()
#
