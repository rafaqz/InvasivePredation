deleteat!(Base.LOAD_PATH, 2:3)

using Revise
using Distributions
using LandscapeChange
using StaticArrays
using Unitful
using Distributions
using ColorSchemes
# using GLMakie
# GLMakie.activate!()
using CairoMakie
CairoMakie.activate!()

using InvasivePredation

basepath = InvasivePredation.basepath

############################################################################
# Calculate mass preferences from papers with both trapped and caught rates

# Data from Childs 1986
(; cat, rodent) = InvasivePredation.load_settings()
(; cat_mass_preference, rodent_stats, rodent_mass_distributions, norway_rat_params, norway_rat_studies) = 
    fit_distributions_to_literature()

# yields = map(rodent_mass_distributions) do rodent
#     (cdf.(rodent, 0:maxmass-1) .- cdf.(rodent, 1:maxmass)) .* 
#     (cdf.(cat_mass_preference, 0:maxmass-1) .- cdf.(cat_mass_preference, 1:maxmass))
# end

# # Makie.lines(cat_mass_preference; strokecolor=:black)
# Makie.lines(yields[1])
# map(Makie.lines!, Base.tail(yields))
# Makie.plot(Normal(1, 1); color=:black)
# map(sum, yields)

fig = let (; norway_rat, black_rat, mouse) = rodent_stats
    # colors = #[:magenta, :cyan, :yellow]
    colors = InvasivePredation.rodent_colors
    fig = Figure(; size=(700, 450));
    kw = (; spinewidth=2, xgridwidth=2, ygridwidth=2)
    ax1 = Axis(fig[1, 1];
        xticks=100:100:700,
        ylabel="Fraction of total",
        kw...
    )
    ax2 = Axis(fig[2, 1]; 
        ylabel="Preferred mass", 
        kw...
    )
    axs = ax1, ax2
    xlims!.(axs, ((0, 700),))
    hidexdecorations!(axs[1])
    hidexdecorations!(axs[2]; ticks=false, ticklabels=false, label=false)
    hideydecorations!.(axs; ticks=false, ticklabels=false, label=false)
    linkxaxes!(axs...)
    alpha = 0.7

    ###################################3
    # A
    nr_colors = ColorSchemes.tempo[[0.15, 0.25]]
    cat_colors = ColorSchemes.solar[[0.75, 0.65]]

    nrs = norway_rat_studies
    gl = nrs.glass.trapped ./ sum(nrs.glass.trapped)
    ch = nrs.childs.trapped ./ sum(nrs.childs.trapped)
    gl_group = [1 for _ in gl]
    ch_group = [2 for _ in ch]
    height = vcat(gl, ch)
    group = vcat(gl_group, ch_group)
    category = vcat(nrs.bin_center_100, nrs.bin_center_100)
    Makie.barplot!(ax1, category, height; 
        label=["Trapped Childs", "Trapped Glass"], 
        color=map(x -> (x, alpha), nr_colors[group]),
        dodge=group,
    )
    Makie.lines!(ax1, nrs.bin_center_100, nrs.glass.killed ./ sum(nrs.glass.killed); 
        label="Predated Glass", 
        color=(cat_colors[1], alpha),
        linewidth=2,
    )
    Makie.lines!(ax1, nrs.bin_center_100, nrs.childs.killed ./ sum(nrs.childs.killed); 
        label="Predated Childs", 
        color=(cat_colors[2], alpha),
        linewidth=2,
    )

    labels = ["Trapped Childs", "Trapped Glass", 
              "Predated Childs", "Predated Glass"]
    elements = vcat(
        [PolyElement(; polycolor=nr_colors[i]) for i in 1:2],
        [LineElement(; linecolor=cat_colors[i]) for i in 1:2],
    )
    Legend(fig[1, 2], elements, labels, ""; 
        framevisible=false,
        halign=:left,
        labelsize=10,
    )

    # B
    Makie.lines!(ax2, cat_mass_preference.distribution; 
        color=InvasivePredation.cat_color, 
        label="Cat prefered mass model",
        linewidth=2,
    )
    Makie.vlines!(ax2, exp(cat_mass_preference.distribution.μ); 
        color=InvasivePredation.cat_color, 
        linestyle=:dash, 
        label="Geometric Mean"
    )

    # Labels
    for (i, text) in enumerate(["A", "B", "C", "D", "E"][1:length(axs)])
        text_ax = Axis(fig[i, 0])
        xlims!(text_ax, (0, 1)); ylims!(text_ax, (0, 1)); simplify!(text_ax)
        text!(text_ax, 0.0, 0.0; text, fontsize=20)
    end
    colsize!(fig.layout, 0, Relative(0.05))

    # Legends
    Legend(fig[2, 2], ax2; 
        framevisible=false, 
        halign=:left,
        labelsize=10,
    )

    fig
end

save(joinpath(basepath, "images/norway_rat_traps_and_kills.png"), fig)
