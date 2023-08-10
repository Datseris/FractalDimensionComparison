export mainplot, mainplot!, estimate_embedding_plot

"""
    rdspl(x, n = 2)
Round `x` to `n` sigdigits for display purposes.
"""
rdspl(x::Real, n = 2) = round(x; digits=n)
rdspl(x::AbstractVector, n = 2) = Tuple((round.(Float64.(x); sigdigits=n)))

function mainplot(args; kwargs...)
    fig, axs = axesgrid(2, 1; sharex = true, kwargs...)
    legendpos = fig[0, :]
    mainplot!(axs, args...; kwargs..., legendpos)
    rowgap!(fig.layout, 10) # spaces the legend table nicely
    return fig
end

function mainplot!(axs, Hs, Cs, eHs, eCs, labels, legendtitle;
        qH = 1, qC = 2, tol = 0.25, tolH = tol, tolC = tol, offsets = zeros(length(Hs)),
        dimension_fit_H = linear_regression_fit_linalg,
        dimension_fit_C = linear_regression_fit_linalg,
        region_choice = :largest, legendpos = nothing, kwargs...
    )

    llines = []
    for j in eachindex(Hs)
        @info labels[j]
        z = offsets[j]
        eH = eHs isa Vector{<:AbstractVector} ? eHs[j] : eHs
        eC = eCs isa Vector{<:AbstractVector} ? eCs[j] : eCs

        # Notice that these plots use Makie's cycling with the current
        # plottheme. This means that colors, markers, and linestyles, are cycled.

        # Entropy
        H = Hs[j]
        x, y = log.(eH), -H
        line = lines!(axs[1], x, y .+ z; alpha = 0.9)
        push!(llines, line)
        # For entropy we always want the largest region because we have not
        # observed any change of slope in neither q≠2 or having noise
        # if region_choice == :largest
            region, d = linear_region(x, y; tol = tolH, warning = false, sat = 0.1)
        # elseif region_choice == :last
        #     lrs, tangents = linear_regions(x, y; tol)
        #     region = lrs[end]
        #     d = tangents[end]
        # end
        Δ, Δ05, Δ95 = dimension_fit_H(x[region], y[region])
        Hlabel = "$(rdspl(Δ05))-$(rdspl(Δ95))"
        scatter!(axs[1], x[[region[1], region[end]]], y[[region[1], region[end]]] .+ z;
            label = Hlabel, alpha = 0.75
        )
        axs[2].xticks = WilkinsonTicks(8; k_min = 6)

        # Correlation sum
        C = Cs[j]
        i = findfirst(c -> c > 0, C)
        if !isnothing(i)
            x, y = log.(eC)[i:end], log.(C)[i:end]
            if region_choice == :largest
                region, d = linear_region(x, y; tol = tolC, warning = false)
            elseif region_choice == :last
                lrs, tangents = linear_regions(x, y; tol = tolC)
                region = lrs[end]
                d = tangents[end]
            end
            Δ, Δ05, Δ95 = dimension_fit_C(x[region], y[region])
            Clabel = "$(rdspl(Δ05))-$(rdspl(Δ95))"
            lines!(axs[2], x, y .+ z; alpha = 0.9)
            scatter!(axs[2], x[[region[1], region[end]]], y[[region[1], region[end]]] .+ z;
                label = Clabel, alpha = 0.75
            )
        end

        axs[2].xticks = WilkinsonTicks(8; k_min = 6)

        yield()
    end

    axs[1].ylabel = L"-H_%$(qH)"
    axs[2].ylabel = L"\log(C_%$(qC))"
    axs[2].xlabel = L"\log(\varepsilon)"

    for ax in axs
        axislegend(ax; labelsize = 20, nbanks = 3, patchsize=(10f0,20), position = :lt, groupgap = 8)
    end

    # Make the informative legend
    if !isnothing(legendpos)
        leg = Legend(legendpos, llines, labels, legendtitle;
            nbanks=3, patchsize=(40f0, 20),
        )
        space_out_legend!(leg)
        space_out_legend!(leg) # not sure why I have to trigger this twice...
    end
    return
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
