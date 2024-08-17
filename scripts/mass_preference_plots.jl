deleteat!(Base.LOAD_PATH, 2:3)

using Revise
using Distributions
using LandscapeChange
using StaticArrays
using Unitful
using GLMakie
using Distributions
# using CairoMakie
# CairoMakie.activate!()

using InvasivePredation

basepath = InvasivePredation.basepath

############################################################################
# Calculate mass preferences from papers with both trapped and caught rates

# Data from Childs 1986


Gnclude("load_data.jl")

(; cat_mass_preference, rodent_stats, params, norway_rat_studies) = 
    optimize_predation_rates_from_literature()

fig = let
    (; norway_rat, black_rat, mouse) = rodent_stats
    fig = Figure(; size=(800, 600));
    ax1 = Axis(fig[1, 1]; ylabel="Fraction trapped")
    ax2 = Axis(fig[2, 1]; ylabel="Fraction of cat kills")
    ax3 = Axis(fig[3, 1]; xlabel="Prey size", ylabel="Reported means")
    ax4 = Axis(fig[4, 1]; ylabel="Probability")
    axs = ax1, ax2, ax3, ax4
    hidexdecorations!.(axs[1:3]; grid=false)
    hidespines!.(axs)
    linkxaxes!(axs...)
    alpha = 0.7
    b1 = Makie.lines!(ax1, norway_rat.bin_center_mass, norway_rat.trap_rate;
        color=(colors[1], alpha), label=rodent_labels[1] * " trapped sizes",
    )
    b2 = Makie.lines!(ax1, black_rat.bin_center_mass, black_rat.trap_rate;
        color=(colors[2], alpha), label=rodent_labels[2] * " trapped sizes",
    )
    b3 = Makie.lines!(ax1, mouse.bin_center_mass, mouse.trap_rate;
        color=(colors[3], alpha), label=rodent_labels[3] * " estimated sizes",
    )
    b4 = Makie.lines!(ax2, params.bin_center_mass, params.killed_rats ./ sum(params.killed_rats);
        color=(colors[1], alpha), label="Norway rat killed sizes",
    )

    d = Makie.density!(ax3, cat_mean_prey_sizes; color=:grey, label="Literature mean prey sizes")
    l = Makie.plot!(ax4, cat_mass_preference; color=:black, label="Preference model")
    axislegend.(axs; position=:rc)
    Makie.xlims!.(axs, ((0, 600),))
    fig
end

save(joinpath(basepath, "images/cat_rodent_predation.png"), fig)

fig = let
    fig = Figure()
    ax = Axis(fig[1, 1])
    nrs = norway_rat_studies
    Makie.barplot!(ax, nrs.bin_center_100, nrs.glass.trapped ./ sum(nrs.glass.trapped); label="Glass trapped")
    Makie.barplot!(ax, nrs.bin_center_100, nrs.childs.trapped ./ sum(nrs.childs.trapped); label="Childs trapped")
    Makie.lines!(ax, nrs.bin_center_100, nrs.glass.killed ./ sum(nrs.glass.killed); label="Glass killed")
    Makie.lines!(ax, nrs.bin_center_25, nrs.childs.killed_25 ./ sum(nrs.childs.killed_25); label="Childs killed")
    axislegend(ax; position=:rt)
    fig
end

save(joinpath(basepath, "images/cat_rodent_predation.png"), fig)
