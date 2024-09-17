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


include(joinpath(basepath, "scripts/load_settings.jl"))

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
    fig = Figure(; size=(800, 600));
    ax1 = Axis(fig[1, 1]; ylabel="Fraction trapped")
    ax2 = Axis(fig[2, 1]; ylabel="LogNormal models")
    ax3 = Axis(fig[3, 1]; ylabel="Fraction of cat kills")
    ax4 = Axis(fig[4, 1]; xlabel="Prey size", ylabel="Reported means")
    ax5 = Axis(fig[5, 1]; ylabel="Probability")
    axs = ax1, ax2, ax3, ax4, ax5
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
    m1 = Makie.lines!(ax2, rodent_mass_distributions.norway_rat;
        color=(colors[1], alpha), label=rodent_labels[1] * " LogNormal model",
    )
    m2 = Makie.lines!(ax2, rodent_mass_distributions.black_rat;
        color=(colors[2], alpha), label=rodent_labels[2] * " LogNormal model",
    )
    m3 = Makie.lines!(ax2, rodent_mass_distributions.mouse;
        color=(colors[3], alpha), label=rodent_labels[3] * " LogNormal model",
    )
    b4 = Makie.lines!(ax3, norway_rat_params.bin_center_mass, norway_rat_params.killed_rats ./ sum(norway_rat_params.killed_rats);
        color=(colors[1], alpha), label="Norway rat killed sizes",
    )

    l = Makie.plot!(ax4, cat_mass_preference; color=:black, label="Cat mass preference model")
    d = Makie.density!(ax5, cat_mean_prey_sizes; color=:grey, label="Literature mean prey sizes")
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

save(joinpath(basepath, "images/norway_rat_traps_and_kills.png"), fig)


