export evtplot
using Statistics

sigr(r, d::Int = 2) = round(r; digits = d)

function evtplot(Ds, labels, legendtitle = "";
        cutoffs = fill(Inf, length(Ds)),
        upperlim = nothing, lowerlim = 0,
        expected = nothing, gap = 0.2, side = :both,
        inner_legend = true,
    )

    fig, axs = axesgrid(1, 1)
    ax = axs[1]
    ax.ylabel = L"\Delta^{(E)}"
    ax.xlabel = "dataset"
    ax.xticks = WilkinsonTicks(length(Ds); k_min = length(Ds))
    violins = []

    for (i, Dloc) in enumerate(Ds)
        # get some statistics
        m = mean(Dloc)
        σ = std(Dloc)
        c = cutoffs[i]
        p = 100count(>(c), Dloc)/length(Dloc)

        lab = "$(sigr(m, 2)) ± $(sigr(σ)) ($(sigr(p, 1))%)"
        horizontal = [i - 0.5 + gap/2, i + 0.5 - gap/2]
        horizontal_small = [i - 0.25 + gap/2, i + 0.25 - gap/2]

        vl = violin!(ax, fill(i, length(Dloc)), Dloc; label = lab, gap, side)
        push!(violins, vl)

        # Line indicating beyond limits
        if c < Inf
            lines!(ax, horizontal_small, [c, c]; color = :red)
        end

        # add mean and expected if there is one
        lines!(ax, horizontal, [m, m]; color = :white, linestyle = :dash)
        if !isnothing(expected)
            e = expected[i]
            lines!(ax, horizontal, [e, e]; color = :white, linestyle = :dot)
        end
    end

    ylims!(ax, (lowerlim, upperlim))
    if inner_legend
        axislegend(ax; labelsize = 20, nbanks = 3, patchsize=(20f0,20), position = :lt, groupgap = 8)
    end

    # Make the informative legend
    leg = Legend(fig[0, :], violins, labels, legendtitle;
        nbanks=3, patchsize=(40f0, 20),
    )
    space_out_legend!(leg)
    space_out_legend!(leg) # not sure why I have to trigger this twice...
    rowgap!(fig.layout, 10)
    display(fig)
    return fig
end