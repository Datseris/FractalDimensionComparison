# Exemplary sets
using DrWatson
@quickactivate :FractalDimensionComparison # re-exports stuff

datas = Vector(undef, 6)
labels = Vector{String}(undef, 6)

N = Int(1e5)
p = 0.99

datas[1] = :roessler_periodic
labels[1] = "periodic"

datas[2] = :henonheiles_quasi
labels[2] = "quasiperiodic"

datas[3] = :koch
labels[3] = "Koch snowflake"

datas[4] = :henon_chaotic
labels[4] = "Hénon"

datas[5] = :lorenz63_chaotic
labels[5] = "Lorenz63"

datas[6] = :kaplanyorke_map
labels[6] = "Kaplan-Yorke"

Dlocs = Vector{Float64}[]

for i in 1:length(datas)
    # These are the parameters that change between the different elements of the plot
    data = datas[i]

    # Here we simply pack all parameters into a dictionary
    # (other parameters are (probably) globals)
    params = @strdict data p N

    if data == :koch
        params["maxk"] = 7
        delete!(params, "N")
    end

    # This is the main call that calculates everything
    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)
end

fig = evtplot(Dlocs, labels, "exemplary sets";
    upperlim = 4, lowerlim = 0.5, cutoffs = [2, 3, 2, 2, 3, 2],
    expected = [1, 2, log(4)/log(3), 1.2, 2.06, 1 - log(2)/log(0.2)]
)

wsave(plotsdir("paper", "evt_analytic"), fig)
