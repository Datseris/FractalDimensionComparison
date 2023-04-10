export mainplot, estimate_embedding_plot

function mainplot(Hs, Cs, eHs, eCs, labels, legendtitle;
        qH = 1, qC = 2, tol = 0.25, offsets = zeros(length(Hs)),
        dimension_fit_H = linear_regression_fit_glm,
        dimension_fit_C = logarithmic_corrected_fit_lsqfit,
    )

    fig, axs = subplotgrid(2, 1; sharex = true)

    llines = []
    for j in 1:length(Hs)
        @info labels[j]
        z = offsets[j]
        eH = eHs isa Vector{<:AbstractVector} ? eHs[j] : eHs
        eC = eCs isa Vector{<:AbstractVector} ? eCs[j] : eCs
        color = Cycled(j)

        C = Cs[j]
        i = findfirst(z -> z > 0, C)
        x, y = log.(eC)[i:end], log.(C)[i:end]
        is, d = linear_region(x, y; tol, warning = false)
        Δ, Δ05, Δ95 = dimension_fit_C(x[is[1]:is[2]], y[is[1]:is[2]])

        # TODO: correct transparency
        Clabel = "$(rdspl(Δ05)), $(rdspl(Δ95))"
        lines!(axs[2], x, y .+ z; color, ls = LINESTYLES[j], alpha = 0.9)
        scatter!(axs[2], x[[is[1], is[end]]], y[[is[1], is[end]]] .+ z; label = Clabel,
            markersize = 10, marker = MARKERS[j], color = color, alpha = 0.75
        )

        H = Hs[j]
        x, y = log.(eH), -H
        line = lines!(axs[1], x, y .+ z; color, ls = LINESTYLES[j], alpha = 0.9)
        push!(llines, line)
        is, d = linear_region(x, y; tol, warning = false)
        Δ, Δ05, Δ95 = dimension_fit_H(x[is[1]:is[2]], y[is[1]:is[2]])
        Hlabel = "$(rdspl(Δ05)), $(rdspl(Δ95))"
        scatter!(axs[1], x[[is[1], is[end]]], y[[is[1], is[end]]] .+ z, label = Hlabel,
        zorder = 2, ms = 10, marker = MARKERS[j], color = color, alpha = 0.75)

        yield()
    end

    axs[1].set_ylabel("\$-H_$(qH)\$")
    axs[2].set_ylabel("\$\\log(C_$(qC))\$")
    axs[2].set_xlabel("\$\\log(\\varepsilon)\$")

    # Make the informative legend
    leg = axs[1].legend(
        handles = llines, labels = labels,
        bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
        ncol=3, mode="expand", borderaxespad=0,
        title = legendtitle
    )

    for ax in axs
        ax.legend(; ncol = 3, handlelength = 0.2)
    end
    axs[1].add_artist(leg)
    fig.tight_layout(pad = 0.25)
    fig.subplots_adjust(top = 0.8)

    return fig, axs
end

function estimate_embedding_plot(x, t = "")
    fig = figure()
    ax = subplot(2,1,1)
    ax.plot(x)
    ax.set_title(t)

    mi = selfmutualinfo(x, 0:100)
    w = estimate_delay(x, "mi_min")
    ax = subplot(2,2,3)
    ax.plot(0:100, mi; label = "mut. inf.")
    ax.axvline(w; color = "C2")
    ax.legend()

    Es = delay_afnn(x, w, 2:10)
    ax = subplot(2,2,4)
    ax.plot(2:10, Es; label = "Cao") # saturates at d=6
    ax.legend()
end
