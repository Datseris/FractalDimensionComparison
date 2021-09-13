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
    X = data_producing_function(; dict2ntuple(params)...)

    method = get(params, "Cmethod", "standard")
    autoexpand = get(params, "autoexpand", false)
    @info "Calculating εs and H..."
    # boxsizes for entropy and corrsum versions
    if haskey(params, "e")
        eH = eC = params["e"]
    else
        # Entropy always gets max range
        eH = estimate_boxsizes(X;
            w = get(params, "w", 1), z = get(params, "z", -2), k = 20,
            autoexpand,
        )
        if method == "bueno"
            P = autoprismdim(X)
            r0 = ChaosTools.estimate_r0_buenoorovio(X, P)
            ε0 = ChaosTools.minimum_pairwise_distance(X)[1]
            if ε0 == 0
                ε0 = r0/(MathConstants.e^3)
            end
            eC = MathConstants.e .^ range(log(ε0)+1, log(r0); length = 12)
        elseif method == "standard"
            eC = eH
        end
    end
    @time H = genentropy.(Ref(X), eH; q = qH)

    # Theiler window
    if !haskey(params, "theiler")
        cols = columns(X)
        theiler = maximum(estimate_delay(x, "mi_min") for x in cols)
        theiler = clamp(theiler, 0, 100) # safety
    else
        theiler = params["theiler"]
    end

    @info "Calculating correlation sum..."
    @time C = boxed_correlationsum(X, eC; q = qC, w = theiler, show_progress=true)
    @pack! params = eH, eC, H, C, theiler
    @info "All done. `make_C_H` is returning now."
    return params
end

export make_C_H
