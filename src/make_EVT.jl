export produce_or_load_EVT

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
