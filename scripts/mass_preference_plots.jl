deleteat!(Base.LOAD_PATH, 2:3)

using Revise
using Distributions
using LandscapeChange
using StaticArrays
using Unitful
using Distributions
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

simplify!(ax; kw...) = (hidedecorations!(ax; kw...); hidespines!(ax))

fig = let (; norway_rat, black_rat, mouse) = rodent_stats
    # colors = #[:magenta, :cyan, :yellow]
    rows = [-1:0, 1:2, 3, 4]
    colors = InvasivePredation.rodent_colors
    fig = Figure(; size=(600, 500));
    kw = (; spinewidth=2, xgridwidth=2, ygridwidth=2)
    ax1 = Axis(fig[rows[1], 1];
        xticks=100:100:700,
        xlabel="Rat mass",
        ylabel="Fraction of\ntotal",
        kw...
    )
    # ax2 = Axis(fig[rows[2], 1]; ylabel="Fraction trapped\nand predated", kw...)
    ax2 = Axis(fig[rows[2], 1]; xlabel="Prey size", ylabel="Probability\ndensity", kw...)
    ax3 = Axis(fig[rows[3], 1]; ylabel="Probability density", kw...)
    axs = ax1, ax2, ax3
    hidexdecorations!.(axs[1:2]; grid=false)
    linkxaxes!(axs...)
    alpha = 0.7
    # A
    nr_colors = ColorSchemes.tempo[[0.2, 0.4]]
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
        color=(ColorSchemes.solar[0.5], alpha),
    )
    Makie.lines!(ax1, nrs.bin_center_25, nrs.childs.killed_25 ./ sum(nrs.childs.killed_25); 
        label="Predated Childs", 
        color=(ColorSchemes.solar[0.8], alpha),
    )
    labels = ["Trapped Childs", "Trapped Glass", 
              "Predated Childs", "Predated Glass"]
    elements = vcat(
        [PolyElement(; polycolor=(nr_colors[i], alpha)) for i in 1:2],
        [LineElement(; linecolor=(nr_colors[i], alpha)) for i in 1:2],
    )

    # B
    # b1, b2, b3 = map(1:3, (norway_rat, black_rat, mouse)) do i, rat
    #     Makie.lines!(ax2, rat.bin_center_mass, rat.trap_rate;
    #         color=(colors[i], alpha), label="Trapped " * rodent.labels[i],
    #     )
    # end
    # m1, m2, m3 = map(1:3, rodent_mass_distributions) do i, rat
    #     Makie.lines!(ax2, rat;
    #         color=(colors[i], alpha), label=rodent.labels[i] * " LogNormal model",
    #     )
    # end
    # b4 = Makie.lines!(ax2, norway_rat_params.bin_center_mass, norway_rat_params.killed_rats ./ sum(norway_rat_params.killed_rats);
    #     color=(ColorSchemes.solar[0.7], alpha), label="Predated Norway rat",
    # )

    # B
    l = Makie.plot!(ax2, cat_mass_preference.distribution; color=:black, label="Prefered mass")
    l = Makie.vlines!(ax2, exp(cat_mass_preference.distribution.μ); color=:grey, linestyle=:dash, label="Mean")
    # C
    d = Makie.density!(ax3, cat.mean_prey_sizes; color=:grey, label="Literature mean\nprey sizes")
    Makie.xlims!.(axs, ((0, 700),))
    # Labels
    for (i, text) in enumerate(["A", "B", "C", "D"][1:length(axs)])
        text_ax = Axis(fig[rows[i], 0])
        xlims!(text_ax, (0, 1)); ylims!(text_ax, (0, 1)); simplify!(text_ax)
        text!(text_ax, 0.0, 0.0; text, fontsize=20)
    end
    colsize!(fig.layout, 0, Relative(0.05))
    # Legends
    legend_titles = ["", "", "Model", "Prey sizes"] 
    for i in 1:length(axs)
        Legend(fig[rows[i], 2], axs[i]; 
            framevisible=false, 
            halign=:left,
            labelsize=10,
            title=legend_titles[i]
        )
    end

    fig
end

save(joinpath(basepath, "images/norway_rat_traps_and_kills.png"), fig)


