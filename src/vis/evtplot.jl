export evtplot
using Statistics

sigr(r, d::Int = 2) = round(r; digits = d)
# sigr(r, d::Int = 2) = @printf("%.$(d)f", a)

function evtplot(Ds, labels, legendtitle = "";
        cutoffs = fill(Inf, length(Ds)),
        upperlim = nothing, lowerlim = 0,
        expected = nothing, gap = 0.2
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

        # I am not sure yet how to make the "beyond limits" results plot for now
        vl = violin!(ax, fill(i, length(Dloc)), Dloc; label = lab, gap,
            # datalimits = (-Inf, c)
        )
        push!(violins, vl)

        # if c < Inf
        #     try
        #         violin!(ax, fill(i, length(Dloc)), Dloc;
        #         datalimits = (c, Inf), color = :transparent,
        #         strokecolor = COLORSCHEME[i], strokewidth = 2,
        #         )
        #     catch
        #     end
        # end

        if c < Inf
            lines!(ax, horizontal_small, [c, c]; color = :red)
        end

        # add mean and expected if there is one
        if !isnothing(expected)
            e = expected[i]
            lines!(ax, horizontal, [m, m]; color = :white, linestyle = :dash)
            lines!(ax, horizontal, [e, e]; color = :white, linestyle = :dot)
        end
    end

    ylims!(ax, (lowerlim, upperlim))
    axislegend(ax; labelsize = 20, nbanks = 3, patchsize=(20f0,20), position = :lt, groupgap = 8)

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