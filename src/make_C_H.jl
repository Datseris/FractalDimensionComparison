export produce_or_load_C_H
using DelayEmbeddings: estimate_delay

function produce_or_load_C_H(params, data; kwargs...)
    # Note: we could make it so there are two different calls to `produce_or_load`
    # so that if H or C change, only that is recomputed. But computing H is so fast,
    # we don't care about that for now...
    output, s = produce_or_load(
        datadir("main"), params, make_C_H;
        filename = params -> savename(params; ignores = ["data"]),
        prefix = string(data), suffix = "jld2", storepatch = false,
        kwargs... # Kwargs are typically `force`
    )
    return output
end


"""
    make_C_H(params)

Generate the values of C (correlation sum) and H (generalized entropy) for the given
parameters. This inclues `qH, qC, N, D` as well as the system name to generate data
from. Used in `produce_or_load` from `DrWatson` and utilizes the `Data` module created
in the source folder.

Important: function deduces suitable boxsizes `e` via the automated pipeline, however
they can also be given as name `"e"` in the `params` dictionary.
By default different boxsizes are calculated for entropy/corrsum, because of the
boxed version requiring smaller boxsizes.
"""
function make_C_H(params)
    @unpack qH, qC, data = params

    # This syntax calls the data generating function and creates the dataset
    # and propagates all keywords of `params` into it (if any can be used, e.g. `D` or `dt`)
    @info "Producing data..."
    data_producing_function = getfield(Data, data)
    @time X = data_producing_function(; dict2ntuple(params)...)

    method = get(params, "Cmethod", "standard")
    autoexpand = get(params, "autoexpand", false)
    # boxsizes for entropy and corrsum versions
    if haskey(params, "e")
        eH = eC = params["e"]
    else
        @info "Calculating boxsizes..."
        # Entropy always gets max range
        @time eH = estimate_boxsizes(X;
            w = get(params, "w", 1), z = get(params, "z", -2), k = 16,
            autoexpand, base = MathConstants.e,
        )
        if method == "bueno"
            P = 2
            r0 = estimate_r0_buenoorovio(X, P)
            ε0 = minimum_pairwise_distance(X)[1]
            if ε0 == 0
                ε0 = r0/(MathConstants.e^3)
            end
            eC = MathConstants.e .^ range(log(ε0)+1, log(r0); length = 12)
        elseif method == "standard"
            eC = eH
        end
    end

    if get(params, "compute_H", true)
        @info "Calculating entropies..."
        # before DynamicalSystems.jl v3.0 this was:
        # @time H = genentropy.(Ref(X), eH; q = qH)
        # now it would be:
        # @time H = [entropy(Renyi(qH, MathConstants.e), ValueHistogram(ε), X) for ε in eH]
        # but we utilize the internal code from FractalDimensions.jl (it's threaded)
        @time H = FractalDimensions._threaded_entropies(X, eH, qH, MathConstants.e, true)
    else
        H = zeros(length(eH))
    end
    # Theiler window
    if !haskey(params, "theiler") || isnothing(params["theiler"])
        cols = eachcol(X)
        theiler = maximum(estimate_delay(x, "mi_min") for x in cols)
        theiler = clamp(theiler, 0, 100) # safety
        @show theiler
    else
        theiler = params["theiler"]
    end

    @info "Calculating correlation sum..."
    use_boxed = get(params, "use_boxed", true)
    P = get(params, "prism", 2)

    if use_boxed
        @time C = boxed_correlationsum(X, eC; P, q = qC, w = theiler, show_progress=true)
    else
        @time C = correlationsum(X, eC; q = qC, w = theiler, show_progress=true)
    end
    @pack! params = eH, eC, H, C, theiler
    @info "All done. `make_C_H` is returning now."
    return params
end

export make_C_H
