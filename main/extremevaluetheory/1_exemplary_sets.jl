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

datas[3] = :lorenz96_chaotic
labels[3] = "1/4 Lorenz96 (D=8)"

datas[4] = :henon_chaotic
labels[4] = "Hénon map"

datas[5] = :towel_chaotic
labels[5] = "Towel map"

datas[6] = :coupled_logistics
labels[6] = "1/4 Coupled logistic (D=8)"

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
    elseif data == :lorenz96_chaotic
        params["D"] = 8
    end

    # This is the main call that calculates everything
    output = produce_or_load_EVT(params, data; force = false)
    @unpack Δloc = output
    push!(Dlocs, Δloc)
end

# Normalize the logistic and lorenz as if they are 2 dimensional
Dlocs[3] = Dlocs[3] ./ 4
Dlocs[6] = Dlocs[6] ./ 4

fig = evtplot(Dlocs, labels, "exemplary sets";
    upperlim = 4, lowerlim = 0.5, cutoffs = [2, 3, 2, 2, 3, 2],
    expected = [1, 2, 6.91/4, 1.26, 2.24, 8/4]
)

display(fig)

wsave(plotsdir("paper", "evt_analytic"), fig)
