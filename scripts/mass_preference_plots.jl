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


include(joinpath(basepath, "scripts/load_settings.jl"))

(; cat_mass_preference, rodent_stats, rodent_mass_distributions, norway_rat_params, norway_rat_studies) = 
    fit_distributions_to_literature()
rodent_stats.black_rat |> pairs
mouse.mean_predation_mass


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
    fig = Figure(; size=(600, 800));
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
    axislegend.(axs; position=:rt)
    Makie.xlims!.(axs, ((0, 600),))
    fig
end

save(joinpath(basepath, "images/cat_rodent_predation.png"), fig)

fig = let
    alpha = 0.6
    fig = Figure(; size=(500, 300))
    ax = Axis(fig[1, 1];
        xticks=100:100:700,
        xlabel="Rat mass",
        ylabel="Fraction of total",
    )
    hidespines!(ax)
    colors = [:violetred1, :deepskyblue]
    nrs = norway_rat_studies
    gl = nrs.glass.trapped ./ sum(nrs.glass.trapped)
    ch = nrs.childs.trapped ./ sum(nrs.childs.trapped)
    gl_group = [1 for _ in gl]
    ch_group = [2 for _ in ch]
    height = vcat(gl, ch)
    group = vcat(gl_group, ch_group)
    category = vcat(nrs.bin_center_100, nrs.bin_center_100)
    Makie.barplot!(ax, category, height; 
        color=map(x -> (x, alpha), colors[group]),
        dodge=group,
    )

    Makie.lines!(ax, nrs.bin_center_100, nrs.glass.killed ./ sum(nrs.glass.killed); 
        label="Glass killed", 
        color=(colors[1], alpha),
    )
    Makie.lines!(ax, nrs.bin_center_25, nrs.childs.killed_25 ./ sum(nrs.childs.killed_25); 
        label="Childs killed", 
        color=(colors[2], alpha),
    )

    labels = ["Trapped Childs 1986", "Trapped Glass 2009", 
              "Killed Childs 1986", "Killed Glass 2009"]
    elements = vcat(
        [PolyElement(; polycolor=(colors[i], alpha)) for i in 1:2],
        [LineElement(; linecolor=(colors[i], alpha)) for i in 1:2],
    )
    title = ""
    Legend(fig[1, 2], elements, labels, title; framevisible=false)

    fig
end

save(joinpath(basepath, "images/norway_rat_traps_and_kills.png"), fig)


