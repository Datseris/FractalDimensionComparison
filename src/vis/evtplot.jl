export evtplot, evtplot!
using Statistics

sigr(r, d::Int = 2) = round(r; digits = d)

function evtplot(args...; kw...)
    fig, axs = axesgrid(1, 1)
    evtplot!(axs[1], args...; kw...)
    return fig
end

function evtplot!(ax::Makie.Axis, Ds, labels, legendtitle = "";
        cutoffs = fill(Inf, length(Ds)),
        upperlim = nothing, lowerlim = 0,
        expected = nothing, gap = 0.2, side = :both,
        inner_legend = true,
        fig = ax.parent, legend_position = fig[0, :],
        xpositions = 1:length(Ds),
        ticknum = 6, obtainer = mean,
    )

    ax.ylabel = L"\Delta^{(E)}"
    ax.xlabel = "dataset"
    ax.xticks = 1:ticknum # WilkinsonTicks(6; k_min = length(Ds))
    violins = []

    for (j, Dloc) in enumerate(Ds)
        i = xpositions[j]
        # get some statistics
        m = mean(Dloc)
        σ = std(Dloc)
        # confidence interval for mean comes from
        # https://www.statology.org/confidence-interval-mean/
        # and same stuff can be found in
        # https://stats.libretexts.org/Courses/Las_Positas_College/Math_40%3A_Statistics_and_Probability/07%3A_Confidence_Intervals_and_Sample_Size/7.02%3A_Confidence_Intervals_for_the_Mean_with_Known_Standard_Deviation\

        # we use the 95% CI
        ci = 1.96*σ/sqrt(length(Dloc))
        left = m - ci
        right = m + ci


        # Various ways to print info about the distributions
        # lab = "$(rdspl(m, 2)) ± $(rdspl(σ)) ($(sigr(p, 1))%)"
        # lab = "($(rdspl(left, 2)), $(rdspl(right, 2))) [$(sigr(p, 1))%]"
        c = cutoffs[j]
        if c < Inf
            p = 100count(>(c), Dloc)/length(Dloc)
            lab = "$(rdspl(m, 2)) [$(round(Int, p))%]"
        else
            lab = "$(rdspl(m, 2))"
        end

        if obtainer == median
            lab = "$(rdspl(median(Dloc), 2))"
        end

        horizontal = [i - 0.5 + gap/2, i + 0.5 - gap/2]
        horizontal_small = [i - 0.25 + gap/2, i + 0.25 - gap/2]

        vl = violin!(ax, fill(i, length(Dloc)), Dloc; label = lab, gap, side)
        push!(violins, vl)

        # Line indicating beyond limits
        if c < Inf
            lines!(ax, horizontal_small, [c, c]; color = :red)
        end

        # add mean and expected if there is one
        lines!(ax, horizontal, [m, m]; color = (:white, 0.75), linestyle = LINESTYLES[5])

        if !isnothing(expected)
            e = expected[j]
            lines!(ax, horizontal, [e, e]; color = (:white, 0.75), linestyle = :dot)
        end
    end

    ylims!(ax, (lowerlim, upperlim))
    if inner_legend
        axislegend(ax; labelsize = 20, nbanks = 3, patchsize=(20f0,20), position = :lt, groupgap = 8)
    end

    # Make the informative legend
    if !isnothing(legend_position)
        leg = Legend(legend_position, violins, labels, legendtitle;
            nbanks=3, patchsize=(40f0, 20),
        )
        space_out_legend!(leg)
        space_out_legend!(leg) # not sure why I have to trigger this twice...
        rowgap!(fig.layout, 10)
    end

end