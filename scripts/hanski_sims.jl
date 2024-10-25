using Revise
using LandscapeChange
using Unitful
using StaticArrays
using StatsBase
using GLMakie
using CairoMakie
using Dispersal
using DynamicGrids
using Setfield
using Rasters
using NCDatasets
using NeutralLandscapes
using Colors
using ColorVectorSpace
using ColorSchemes

using InvasivePredation

basepath = InvasivePredation.basepath
alphapath = joinpath(basepath, "tables/alpha")
(; cat, rodent) = s = InvasivePredation.load_settings()
(; cat_mass_preference, rodent_stats, norway_rat_params, norway_rat_studies) =
    fit_distributions_to_literature()
(; stochastic_rates, max_yield_fraction) = InvasivePredation.get_max_yield_fraction()
(; assimilated_energy, individuals_per_cat) = InvasivePredation.get_cat_energetics(cat, rodent, rodent_stats)

t = u"yr" / 12
predation_rates = NamedVector(map(s -> s.predation_rate, rodent_stats))
# Does this make sense ? What is kg and km here exactly
Km_pred = 3.4 * 10^1 * u"kg*km^-2" # (from Mulder 2016)
Ds = uconvert.(u"ha^-1", predation_rates.norway_rat ./ predation_rates * Km_pred / rodent_stats.norway_rat.mean_mass)
cs = individuals_per_cat # Base intake reequirement
ys = max_yield_fraction ./ 2 # Prey yields
rs = rodent.rmax
v = cat.rmax
e = cat.energy_intake
Es = assimilated_energy
d_high = 0.2 / t
# Ncrit = 4 / u"d"

# Define α parameters based on energy intake
# Using mouse scaling for all currently, rat numbers were very wrong
metabolic_rate = (;
    norway_rat=637.0u"kJ*kg^-(3/4)*d^-1",
    black_rat=637.0u"kJ*kg^-(3/4)*d^-1",
    mouse=637.0u"kJ*kg^-(3/4)*d^-1",
)


# Parameter estimates

lc_energy_available = NamedVector(;
    forest=5000u"kJ*d^-1",
    cleared=5000u"kJ*d^-1",
    urban=15000u"kJ*d^-1",
)

lc_exploitation_capacity = NamedVector(;
    forest=NamedVector(;
        norway_rat = 0.5,
        black_rat  = 1.0, # Aboreal
        mouse      = 0.7,
    ),
    cleared=NamedVector(;
        norway_rat = 0.7, # Avoid open fields
        black_rat  = 0.7, # Avoid open fields
        mouse      = 1.0,
    ),
    urban=NamedVector(;
        norway_rat = 1.0,
        black_rat  = 1.0,
        mouse      = 1.0,
    ),
)

# Diet overlap 
nr_br = 0.6
nr_m = 0.4
br_m = 0.4

diet_overlap = NamedVector(;
    norway_rat = (black_rat=nr_br, mouse=nr_m),
    black_rat  = (norway_rat=nr_br, mouse=br_m),
    mouse      = (norway_rat=nr_m, black_rat=br_m),
)

aggression_exponent = 1//3

# Calculated
energy_intake = map(metabolic_rate, rodent_stats) do mr, rs
    uconvert(u"kJ*d^-1", mr * uconvert(u"kg", rs.mean_mass)^(3/4))
end |> NamedVector

rodent_keys = NamedVector{propertynames(rodent_masses)}(propertynames(rodent_masses))
competitor_names = NamedVector(
    norway_rat=(:black_rat, :mouse),
    black_rat=(:norway_rat, :mouse),
    mouse=(:norway_rat, :black_rat)
)
rodent_masses = NamedVector(map(r -> r.mean_mass, rodent_stats))
aggression = (rodent_masses ./ rodent_masses.mouse) .^ aggression_exponent
α = map(diet_overlap, rodent_keys, energy_intake, aggression, competitor_names) do dio, k, ei, c, others
    map(others) do o
        energy_intake[o] / ei  * dio[o] .* max(1.0, aggression[o] / aggression[k])
    end
end; pairs(α)
lc_ks = map(lc_energy_available, lc_exploitation_capacity) do ea, exploitation_capacity
    map(exploitation_capacity, energy_intake) do ec, ei
        upreferred(ec .* ea ./ ei) / u"ha"
    end
end
lc_labels = ["Forest", "Cleared", "Urban"]
lc_names = NamedVector{propertynames(lc_ks)}(propertynames(lc_ks))
scenario_labels = 'A':'E'

##############################################################################3
# Spatially-Implicit model

# Init conditions
P0 = 0.01 * u"ha^-1"
ks = lc_ks.forest
Ns0 = ks ./ 2
tspan = 1u"d":1u"d":convert(typeof(1.0u"d"), 1u"yr")
stochasticity = Param(0.05; bounds=(0.0, 0.2), label="stochasticity")

model = (; P0, Ns0, ks, Ds, Es, ys, α, cs, rs, d_high, v, e, t, tspan, stochasticity, cat)

function plot_populations!(fig, model, lc_ks, i; 
    lc=i, 
    hideticks=true, 
    label=false,
)
    cat_color = InvasivePredation.cat_color
    rodent_colors = InvasivePredation.rodent_colors
    m = @set model.ks = lc_ks[lc]

    l = ('A':'E')[i] * ": " * lc_labels[lc]
    label && label!(fig, l, 0, i)
    # Populations
    ax = Axis(fig[1, i]; 
        yscale=log10, 
        xticks=[100, 200, 300],
        yticks=[0.01, 0.1, 1, 10, 100],
    )
    ylims!(ax, (0.001, 500))
    results = hanski_sim(m)
    log_fudge = (0.000000001 * oneunit(last(last(last(results)))))
    map(3:-1:1) do r
        rodent_timeline = getindex.(last(results), r) .+ log_fudge
        Makie.lines!(ax, ustrip.(u"d", tspan), rodent_timeline; 
            label=rodent.labels[r], 
            color=rodent_colors[r],
            linewidth=2,
        )
    end
    cat_timeline = first(results) .+ log_fudge
    Makie.lines!(ax, ustrip.(u"d", tspan), cat_timeline; 
        label="Cat", 
        color=cat_color,
        linewidth=2,
    )
    hideydecorations!(ax; ticks=hideticks, ticklabels=hideticks)
    hidexdecorations!(ax; label=false, ticks=false, ticklabels=false)
    hidespines!(ax)
    return ax
end

CairoMakie.activate!()
fig = let
    fig = Figure(; size=(800, 250));
    axs = map(1:3) do i
        plot_populations!(fig, model, lc_ks, i; hideticks=i != 1, label=true)
    end
    Legend(fig[1, 4], axs[1]; framevisible=false)
    rowgap!(fig.layout, 1, 0)
    rowsize!(fig.layout, 0, Relative(0.1))
    fig
end

save(joinpath(basepath, "images/si_sim.png"), fig)
save(joinpath(basepath, "images/si_sim.svg"), fig)

fig = let
    fig = Figure(; size=(800, 250));
    ax = plot_populations!(fig, model, lc_ks, 1; hideticks=true, lc=3)
    Legend(fig[1, 2], ax; framevisible=false)
    fig
end

save(joinpath(basepath, "images/si_sim_mini.png"), fig)
save(joinpath(basepath, "images/si_sim_mini.svg"), fig)

# Live interaction to find α parameters
# mm = MakieModel(model) do layout, m
#     ax1 = Axis(layout[1, 1])
#     ax2 = Axis(layout[2, 1])
#     res = lift(hanski_sim, m)
#     P_timeline = lift(first, res)
#     Rn_timeline = lift(r -> getindex.(last(r), 1), res)
#     Rr_timeline = lift(r -> getindex.(last(r), 2), res)
#     Mm_timeline = lift(r -> getindex.(last(r), 3), res)
#     Makie.lines!(ax1, Rn_timeline; label=rodent.labels[1])
#     Makie.lines!(ax1, Rr_timeline; label=rodent.labels[2])
#     Makie.lines!(ax1, Mm_timeline; label=rodent.labels[3])
#     Makie.lines!(ax2, P_timeline; label="Cat")
#     linkxaxes!(ax1, ax2)
#     axislegend(ax1)
#     axislegend(ax2)
# end


##############################################################################3
# Spatially-Explicit model

hanski_rule = let lc=Aux{:lc}(), lc_ks=stripparams(lc_ks), model=model
    Cell{Tuple{:cats,:rodents}}() do data, (P, Ns), I
        lc_i = get(data, lc, I)
        lc_i == 0 && return zero(P), zero(Ns)
        ks = lc_ks[lc_i]
        Ns1 = hanski_prey_timestep(P, Ns, ks, model.α, model)::typeof(Ns)
        P1 = hanski_predator_timestep(P, Ns, model)::typeof(P)
        return P1, Ns1
    end
end

cat_spread_rule = OutwardsDispersal{:cats}(;
    stencil=Moore(3),
    formulation=ExponentialKernel(Param(0.2, label="cat λ", bounds=(0.00001, 10.0))),
    maskbehavior=Dispersal.CheckMaskEdges()
)

rodent_λ = NamedVector{propertynames(Ns0)}((
    Param(0.2, label="R.n λ", bounds=(0.00001, 1.0)),
    Param(0.2, label="R.r λ", bounds=(0.00001, 1.0)),
    Param(0.2, label="M.n λ", bounds=(0.00001, 1.0)),
))
ExponentialKernel([1.0, 2.0])(10.0)

rodent_spread_rule = InwardsDispersal{:rodents}(;
    stencil=Moore(3),
    formulation=ExponentialKernel(rodent_λ),
    # maskbehavior=Dispersal.CheckMaskEdges()
)

# rodent_allee_rule = AlleeExtinction{:rodents}(
#     NamedVector(; norway_rat=0.4, black_rat=0.4, mouse=0.4) .* u"ha^-1"
# )

ruleset = Ruleset(hanski_rule, rodent_spread_rule, cat_spread_rule)#, rodent_allee_rule)

rodents = Raster(fill(lc_ks.forest ./ 2, X(64), Y(64)))
cats = map(_ -> 0.01oneunit(lc_ks.forest[1]), rodents)
tspan = 1:1000
mpd = Raster(rand(MidpointDisplacement(0.35), dims(rodents)))
a = quantile(vec(mpd), 1/4)
b = quantile(vec(mpd), 3/4)
classes = (<(a) => 1, a .. b => 2, >(b) => 3)
lc = Rasters.classify(mpd, classes; others=0, missingval=0)
init = (; cats, rodents)

# Basic test run
output = ResultOutput(init; tspan=1:10, aux=(; lc=parent(lc)), boundary=Wrap())
sim!(output, ruleset; printframe=true)

#####################################
# Live simulation

# GLMakie.activate!()
# output = MakieOutput(init;
#     ruleset, tspan, fps=100,
#     printframe=true,
#     aux=(; lc=parent(lcs[lc_mode])),
#     store=false,
#     ncolumns=2,
# ) do x
#     colors = (RGBA(0.5, 0.5, 0.0), RGBA(0.5, 0.0, 0.5), RGBA(0.0, 0.5, 0.5))
#     inds = [(1, 1), (2, 1), (2, 2), (2, 3)]
#     I = inds[1]

#     # Landcover
#     lc_ax = Axis(x.layout[1, 1]; xlabel="land cover")
#     image!(lc_ax, lcs[lc_mode]; colormap=:amp, interpolate=false, colorrange=(0, 4))
#     # Cats
#     cat_ax = Axis(x.layout[1, 2]; xlabel="cat")
#     c = image!(cat_ax, x.frame.cats; colormap=:amp, interpolate=false, colorrange=(0.0, 0.07))
#     Colorbar(x.layout[1, 2, Right()], c; flipaxis=false)
#     # Rodents
#     rodent_axs = map(1:length(lc_ks[1]), propertynames(lc_ks[1]), maximum(lc_ks)) do i, name, k
#         I = inds[i + 1]
#         rodent_ax = Axis(x.layout[I...]; xlabel=string(name))
#         A = rebuild(getindex.(rodents, i); missingval=0.0)
#         rodent_obs = Observable{Raster}(A)
#         on(x.frame.rodents) do rodents
#             rodent_obs[] .= getindex.(rodents, i)
#             notify(rodent_obs)
#         end
#         r = image!(rodent_ax, rodent_obs; colormap=:tempo, interpolate=false, colorrange=(0.0, k))
#         Colorbar(x.layout[I..., Right()], r; flipaxis=false)
#         rodent_ax
#     end
#     linkaxes!(cat_ax, lc_ax, rodent_axs...)
#     nothing
# end

#######################################
# Figures
CairoMakie.activate!()

# Run simulations
outputs = let tspan=1:500, aux=(; lc=parent(lc))
    init1 = (; cats, rodents)
    init2 = (; cats, rodents=map(x -> x .* NamedVector(norway_rat=0, black_rat=1, mouse=1), rodents))
    init3 = (; cats, rodents=map(x -> x .* NamedVector(norway_rat=1, black_rat=0, mouse=1), rodents))
    init4 = (; cats, rodents=map(x -> x .* NamedVector(norway_rat=1, black_rat=1, mouse=0), rodents))
    init5 = (; cats=rebuild(cats, parent(cats .* 0)), rodents)
    outputs = (
        ResultOutput(init1; tspan, aux),
        ResultOutput(init2; tspan, aux),
        ResultOutput(init3; tspan, aux),
        ResultOutput(init4; tspan, aux),
        ResultOutput(init5; tspan, aux),
    )
    foreach(outputs) do output
        sim!(output, ruleset; printframe=true, boundary=Wrap())
    end
    outputs
end

fig = let ylabelsize=20
    fig = Figure(; size=(1250, 1200))
    grays = cgrad(ColorSchemes.grayC[((0:1:3) ./ 3)[1:3]], 3; categorical=true)
    r_maxs = [200.0, 200.0, 200.0]
    r_ticks = [0:50:150, 0:50:150, 0:50:150]
    colorbar_kw = (;
        flipaxis=true,
        spinewidth=0,
        width=20,
    )
    all_axs = map(enumerate(outputs)) do (j, output)
        # Labels
        text_ax = Axis(fig[0, j]; axis=1)
        xlims!(text_ax, (0, 1))
        ylims!(text_ax, (0, 1))
        simplify!(text_ax)
        text!(text_ax, 0.0, 0.0; text=scenario_labels[j], fontsize=20)
        # Landcover
        lc_ax = Axis(fig[1, j]; ylabel="Land Cover", ylabelsize, aspect=1)
        lc_im = image!(lc_ax, lc; colormap=grays, interpolate=false, colorrange=(0.5, 3.5))
        j == length(outputs) && Colorbar(fig[1, j + 1], lc_im;
            ticks=(1:1:3, lc_labels),
            colorbar_kw...
        )
        # Rodents
        rodent_axs = map(1:length(lc_ks[1]), propertynames(lc_ks[1]), r_maxs) do i, name, max
            rodent_ax = Axis(fig[i + 1, j];
                ylabel=titlecase(replace(string(name), "_"=> " ")),
                ylabelsize, 
                aspect=1
            )
            rA = output[end].rodents
            r = image!(rodent_ax, getindex.(rA, i);
                colormap=:tempo,
                interpolate=false,
                colorrange=(0.0, max),
            )
            j == length(outputs) && Colorbar(fig[i + 1, j + 1], r;
                ticks=r_ticks[i],
                colorbar_kw...
            )
            rodent_ax
        end
        # Cats
        cat_ax = Axis(fig[5, j]; ylabel="Cat", ylabelsize, aspect=1)
        cA = output[end].cats
        c = image!(cat_ax, cA;
            colormap=Reverse(:solar),
            interpolate=false,
            colorrange=(0.0, 0.07),
        )
        j == length(outputs) && Colorbar(fig[5, j + 1], c;
            ticks=collect(0.01:0.01:0.06),
            colorbar_kw...
        )
        axs = (cat_ax, lc_ax, rodent_axs...)
        # No boxes or ticks
        j == 1 ? simplify!.(axs; label=false) : simplify!.(axs)
        axs
    end |> Iterators.flatten |> collect
    # Zoom everything together
    linkaxes!(all_axs...)
    rowsize!(fig.layout, 0, Relative(0.03))
    fig
end

save(joinpath(basepath, "images/se_sim.png"), fig)
save(joinpath(basepath, "images/se_sim.svg"), fig)

fig = let ylabelsize=20
    fig = Figure(; size=(1000, 220))
    grays = cgrad(ColorSchemes.grayC[((0:1:3) ./ 3)[1:3]], 3; categorical=true)
    output = outputs[1]
    # Landcover
    lc_ax = Axis(fig[1, 1]; aspect=1)
    lc_im = image!(lc_ax, lc; colormap=grays, interpolate=false, colorrange=(0.5, 3.5))
    # Rodents
    rodent_axs = map(1:3) do i
        ax = Axis(fig[1, 1 + i]; aspect=1)
        A = rebuild(getindex.(output[end].rodents, i); missingval=0.0u"ha^-1")
        r = image!(ax, A;
            colormap=:tempo,
            interpolate=false,
            colorrange=(0.0, 150.0),
        )
        ax
    end
    # Cats
    cat_ax = Axis(fig[1, 2 + length(rodent_axs)]; aspect=1)
    c = image!(cat_ax, rebuild(output[end].cats; missingval=0.0u"ha^-1");
        colormap=Reverse(:solar),
        interpolate=false,
        colorrange=(0.0, 0.07),
    )
    axs = (cat_ax, lc_ax, rodent_axs...)
    simplify!.(axs)
    fig
end

save(joinpath(basepath, "images/se_sim_mini.png"), fig)
save(joinpath(basepath, "images/se_sim_mini.svg"), fig)

# fig

