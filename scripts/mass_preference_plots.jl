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

fig = let (; norway_rat, black_rat, mouse) = rodent_stats
    # colors = #[:magenta, :cyan, :yellow]
    colors = [:red, :lightblue, :yellow]
    fig = Figure(; size=(600, 600));
    ax1 = Axis(fig[1:2, 1]; ylabel="Fraction trapped\nand predated")
    ax2 = Axis(fig[3, 1]; xlabel="Prey size", ylabel="Probability\ndensity")
    ax3 = Axis(fig[4, 1]; ylabel="Probability\ndensity")
    ax4 = Axis(fig[-1:0, 1];
        xticks=100:100:700,
        xlabel="Rat mass",
        ylabel="Fraction of\ntotal",
    )
    axs = ax1, ax2, ax3, ax4
    hidexdecorations!.(axs[1:2]; grid=false)
    hidespines!.(axs)
    linkxaxes!(axs...)
    alpha = 0.7
    b1, b2, b3 = map(1:3, (norway_rat, black_rat, mouse)) do i, rat
        Makie.lines!(ax1, rat.bin_center_mass, rat.trap_rate;
            color=(colors[i], alpha), label=rodent.labels[i] * " trapped sizes",
        )
    end
    # m1, m2, m3 = map(1:3, rodent_mass_distributions) do i, rat
    #     Makie.lines!(ax2, rat;
    #         color=(colors[i], alpha), label=rodent.labels[i] * " LogNormal model",
    #     )
    # end
    b4 = Makie.lines!(ax1, norway_rat_params.bin_center_mass, norway_rat_params.killed_rats ./ sum(norway_rat_params.killed_rats);
        color=(:black, alpha), label="Norway rat predated mass",
    )

    l = Makie.plot!(ax2, cat_mass_preference.distribution; color=:black, label="Model of prefered rodent mass")
    l = Makie.vlines!(ax2, exp(cat_mass_preference.distribution.μ); color=:grey, linestyle=:dash, label="Mean preferred mass")
    d = Makie.density!(ax3, cat.mean_prey_sizes; color=:grey, label="Literature mean prey sizes")
    Makie.xlims!.(axs, ((0, 700),))

    alpha = 0.6
    hidespines!(ax4)
    colors = [:darkred, :red]
    nrs = norway_rat_studies
    gl = nrs.glass.trapped ./ sum(nrs.glass.trapped)
    ch = nrs.childs.trapped ./ sum(nrs.childs.trapped)
    gl_group = [1 for _ in gl]
    ch_group = [2 for _ in ch]
    height = vcat(gl, ch)
    group = vcat(gl_group, ch_group)
    category = vcat(nrs.bin_center_100, nrs.bin_center_100)
    Makie.barplot!(ax4, category, height; 
        label=["Trapped Childs", "Trapped Glass"], 
        color=map(x -> (x, alpha), colors[group]),
        dodge=group,
    )

    Makie.lines!(ax4, nrs.bin_center_100, nrs.glass.killed ./ sum(nrs.glass.killed); 
        label="Predated Glass", 
        color=(colors[1], alpha),
    )
    Makie.lines!(ax4, nrs.bin_center_25, nrs.childs.killed_25 ./ sum(nrs.childs.killed_25); 
        label="Predated Childs", 
        color=(colors[2], alpha),
    )

    labels = ["Trapped Childs 1986", "Trapped Glass 2009", 
              "Predated Childs 1986", "Predated Glass 2009"]
    elements = vcat(
        [PolyElement(; polycolor=(colors[i], alpha)) for i in 1:2],
        [LineElement(; linecolor=(colors[i], alpha)) for i in 1:2],
    )
    title = ""
    # Legend(fig[-1:0, 1], elements, labels, title; framevisible=false)
    axislegend.(axs[1:4]; position=:rt)

    fig
end

save(joinpath(basepath, "images/norway_rat_traps_and_kills.png"), fig)


