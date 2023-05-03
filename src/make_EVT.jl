export produce_or_load_EVT, produce_or_load_pointwise

function produce_or_load_EVT(params, data; kwargs...)

    output, s = produce_or_load(
        datadir("evt"), params, make_EVT;
        filename = params -> savename(params; ignores = ["data"]),
        prefix = "evt_"*string(data), suffix = "jld2", storepatch = false,
        kwargs... # Kwargs are typically `force`
    )
    return output
end

function make_EVT(params)
    @unpack p, data = params

    @info "Producing data..."
    data_producing_function = getfield(Data, data)
    X = data_producing_function(; dict2ntuple(params)...)
    estimator = get(params, "estimator", :mm)
    allocate_matrix = length(X) ≤ 10_000

    @info "Calculating EVT local dims..."
    Δloc, = extremevaltheory_dims_persistences(X, p;
        allocate_matrix, estimator, compute_persistence = false
    )
    @pack! params = Δloc
    return params
end

function produce_or_load_pointwise(params, data; kwargs...)

    function make_pairwise(params)
        @info "Producing data..."
        data_producing_function = getfield(Data, data)
        X = data_producing_function(; dict2ntuple(params)...)

        es = estimate_boxsizes(X)
        @info "Calculating pointwise local dims..."
        Cs = pointwise_correlationsums(X, es)

        x = log.(es)
        function local_slope(x, C)
            L = length(x)
            i = findfirst(c -> c > 0, C)
            z = view(x, i:L)
            y = log.(C)[i:end]
            return linear_region(z, y; warning = false)[2]
        end
        Dlocs = map(C -> local_slope(x, C), Cs)

        params["Cs"] = Cs
        params["es"] = es
        params["Ds"] = Dlocs
        return params
    end

    output, s = produce_or_load(
        datadir("evt"), params, make_pairwise;
        filename = params -> savename(params; ignores = ["data"]),
        prefix = "pointwise_"*string(data), suffix = "jld2", storepatch = false,
        kwargs... # Kwargs are typically `force`
    )
    return output
end