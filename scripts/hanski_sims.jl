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

# aggfactor = 4
# lc_path = joinpath(basepath, "data/lc_predictions_mus.nc")
# lc = RasterStack(lc_path) |>
#     x -> Rasters.maybeshiftlocus(Rasters.Start(), x) |>
#     x -> Rasters.set(x, Ti => Int.(lookup(x, Ti))) |>
#     x -> Rasters.aggregate(Rasters.Start(), x, (X(aggfactor), Y(aggfactor))) |>
#     x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |> NamedTuple
# lc = lc.native .* 1 .+ lc.cleared .* 2 .+ lc.abandoned .* 3  .+ lc.forestry .* 3 .+ lc.urban .* 4
# m = lc .!= 0

t = u"yr" / 12
predation_rates = NamedVector(map(s -> s.predation_rate, rodent_stats))
# Hanski D parameters taken from cat mass preference model
# These are relative, absolute value doesn't matter
Ds = predation_rates ./ predation_rates[1]
cs = individuals_per_cat # Base intake reequirement
ys = max_yield_fraction ./ 2 # Prey yields
v = cat.rmax * u"yr^-1" # Predator rmax
e = cat.energy_intake
Es = assimilated_energy
rs = rodent.rmax
d_high = 0.2 / t
# Ncrit = 4 / u"d"


# Define α parameters based on energy intake and diet overlap
metabolic_rate = (;
    norway_rat=211.0u"kJ*kg^-(3/4)*d^-1",
    black_rat=211.0u"kJ*kg^-(3/4)*d^-1",
    mouse=211.0u"kJ*kg^-(3/4)*d^-1",
)

energy_intake = map(metabolic_rate, rodent_stats) do mr, rs
    uconvert(u"kJ/d", mr * (rs.mean_mass * u"g"^(3/4)))
end |> NamedVector

lc_energy_available = NamedVector(;
    native=5000u"kJ*d^-1",
    cleared=5000u"kJ*d^-1",
    urban=10000u"kJ*d^-1",
)

lc_exploitation_capacity = NamedVector(;
    native=NamedVector(;
        norway_rat = 0.4,
        black_rat  = 0.8,
        mouse      = 0.7,
    ),
    cleared=NamedVector(;
        norway_rat = 0.5,
        black_rat  = 0.5,
        mouse      = 0.8,
    ),
    urban=NamedVector(;
        norway_rat = 0.8,
        black_rat  = 0.6,
        mouse      = 0.8,
    ),
)

lc_energy_available = NamedVector(;
    native=8000u"kJ*d^-1",
    cleared=8000u"kJ*d^-1",
    urban=24000u"kJ*d^-1",
)

#=
native: 8000 kJ*d^-1, cleared: 8000 kJ*d^-1, urban: 40000 kJ*d^-1.
native: (0.3, 0.8, 0.6), cleared: (0.5, 0.5, 0.5), urban: (0.8, 0.6, 0.6,)
=#
nr_br = 0.7
nr_m = 0.6
br_m = 0.6

diet_overlap = NamedVector(;
    norway_rat = (nr_br, nr_m),
    black_rat  = (nr_br, br_m),
    mouse      = (nr_m, br_m),
)

competitors = NamedVector(norway_rat=(2, 3), black_rat=(1, 3), mouse=(1, 2))
α = map(diet_overlap, energy_intake, competitors) do di_o, ei, (a, b)
    (energy_intake[a] / ei * di_o[1], energy_intake[b] / ei  * di_o[2])
end

lc_ks = map(lc_energy_available, lc_exploitation_capacity) do ea, exploitation_capacity
    map(exploitation_capacity, energy_intake) do ec, ei
        upreferred(ec .* ea ./ ei)
    end
end
lc_labels = ["Native", "Cleared", "Urban"]
scenario_labels = ["A", "B", "C"]

##############################################################################3
# Spatially-Implicit model

# Init conditions
P = 0.0001
Ns = rodent.carrycap .* u"ha" ./ 2
tspan = 1u"d":1u"d":convert(typeof(1.0u"d"), 3u"yr")
stochasticity = Param(0.05; bounds=(0.0, 0.2), label="stochasticity")
lc_names = NamedVector{keys(lc_ks)}(keys(lc_ks))
ks = lc_ks.native

model = (; P, Ns, ks, Ds, Es, ys, α, cs, rs, d_high, v, e, t, tspan, stochasticity, cat)

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

simplify!(ax; kw...) = (hidedecorations!(ax; kw...); hidespines!(ax))

CairoMakie.activate!()
fig = let
    cat_color = InvasivePredation.cat_color
    rodent_colors = InvasivePredation.rodent_colors
    fig = Figure(; size=(800, 250));
    for i in 1:3
        m = @set model.ks = lc_ks[i]
        # Label
        text_ax = Axis(fig[0, i])
        xlims!(text_ax, (0, 1))
        ylims!(text_ax, (0, 1))
        simplify!(text_ax)
        text!(text_ax, 0.0, 0.0; text=scenario_labels[i], fontsize=20)
        # Populations
        ax = Axis(fig[1, i]; yscale=log10, xlabel=lc_labels[i])
        results = hanski_sim(m)
        map(3:-1:1) do r
            rodent_timeline = getindex.(last(results), r) .+ 0.000000001
            Makie.lines!(ax, tspan, rodent_timeline; label=rodent.labels[r], color=rodent_colors[r])
        end
        cat_timeline = first(results) .+ 0.000000001
        Makie.lines!(ax, tspan, cat_timeline; label="Cat", color=cat_color)
        hideydecorations!(ax; ticks=i != 1, ticklabels=i != 1) 
        hidexdecorations!(ax; label=false, ticks=false, ticklabels=false) 
        hidespines!(ax)
        if i == 3 
            Legend(fig[1, 4], ax; framevisible=false)
        end
    end
    rowsize!(fig.layout, 0, Relative(0.2))
    fig
end


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

rodent_λ = NamedVector{propertynames(Ns)}((
    Param(0.5, label="R.n λ", bounds=(0.00001, 1.0)),
    Param(0.3, label="R.r λ", bounds=(0.00001, 1.0)),
    Param(0.2, label="M.n λ", bounds=(0.00001, 1.0)),
))
ExponentialKernel([1.0, 2.0])(10.0)

rodent_spread_rule = InwardsDispersal{:rodents}(;
    stencil=Moore(3),
    formulation=ExponentialKernel(rodent_λ),
    # maskbehavior=Dispersal.CheckMaskEdges()
)

rodent_allee_rule = AlleeExtinction{:rodents}(NamedVector(; norway_rat=0.1, black_rat=0.4, mouse=0.5))

ruleset = Ruleset(hanski_rule, rodent_spread_rule, cat_spread_rule, rodent_allee_rule)

rodents = Raster(fill(lc_ks.native ./ 2, X(150), Y(300)))
cats = map(_ -> 0.01, rodents)
tspan = 1:1000
quantiles = (
    (0.60, 0.80),
    (0.20, 0.80),
    (0.20, 0.40),
)
lcs = map(quantiles) do (qa, qb)
    mpd = Raster(rand(MidpointDisplacement(0.35), dims(rodents)))
    a = quantile(vec(mpd), qa)
    b = quantile(vec(mpd), qb)
    classes = (<(a) => 1, a .. b => 2, >(b) => 3)
    Rasters.classify(mpd, classes; others=0, missingval=0)
end
tspan = 1:1000
init = (; cats, rodents)
lc_mode = 3

# Basic test run
output = ResultOutput(init; tspan=1:10, aux=(; lc=parent(lcs[lc_mode])), boundary=Wrap())
sim!(output, ruleset; printframe=true)

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

# Figure

CairoMakie.activate!()

# Run simulations
outputs = map(1:3) do i
    output = ResultOutput(init; tspan=1:500, aux=(; lc=parent(lcs[i])), boundary=Wrap())
    sim!(output, ruleset; printframe=true)
end

fig = let ylabelsize=20
    fig = Figure(; size=(1000, 1200))
    # maxs = mapreduce((a, x) -> map(max, a, x), outputs) do output
    #     reduce(output[end].rodents) do acc, x
    #         map(max, acc, x)
    #     end
    # end
    grays = cgrad(ColorSchemes.grayC[((0:1:3) ./ 3)[1:3]], 3; categorical=true)
    r_maxs = [60.0, 60.0, 300.0]
    r_ticks = [0:10:50, 0:10:50, 0:50:250]
    all_axs = map(enumerate(outputs)) do (j, output)
        # Labels
        text_ax = Axis(fig[0, j])
        xlims!(text_ax, (0, 1))
        ylims!(text_ax, (0, 1))
        simplify!(text_ax)
        text!(text_ax, 0.0, 0.0; text=labels[j], fontsize=20)
        # Landcover
        lc_ax = Axis(fig[1, j]; ylabel="Land Cover", ylabelsize)
        lc_im = image!(lc_ax, lcs[j]; colormap=grays, interpolate=false, colorrange=(0.5, 3.5))
        j == 3 && Colorbar(fig[1, j, Right()], lc_im;
            flipaxis=true,
            spinewidth=0,
            ticks=(1:1:3, lc_labels),
        )
        # Rodents
        rodent_axs = map(1:length(lc_ks[1]), propertynames(lc_ks[1]), r_maxs) do i, name, max
            I = 1 + i, j
            rodent_ax = Axis(fig[I...];
                ylabel=titlecase(replace(string(name), "_"=> " ")),
                ylabelsize,
            )
            A = rebuild(getindex.(output[end].rodents, i); missingval=0.0)
            r = image!(rodent_ax, A;
                colormap=:tempo,
                interpolate=false,
                colorrange=(0.0, max),
            )
            j == 3 && Colorbar(fig[I..., Right()], r;
                flipaxis=true,
                ticks=r_ticks[i],
                spinewidth=0,
            )
            rodent_ax
        end
        # Cats
        cat_ax = Axis(fig[5, j]; ylabel="Cat", ylabelsize)
        c = image!(cat_ax, output[end].cats;
            colormap=Reverse(:solar),
            interpolate=false,
            colorrange=(0.0, 0.07),
        )
        j == 3 && Colorbar(fig[5, j, Right()], c;
            flipaxis=true,
            spinewidth=0,
            ticks=collect(0.01:0.01:0.06),
        )
        axs = (cat_ax, lc_ax, rodent_axs...)
        # No boxes or ticks
        j == 1 ? simplify!.(axs; label=false) : simplify!.(axs)
        axs
    end |> Iterators.flatten |> collect

    rowsize!(fig.layout, 0, Relative(0.05))
    ax = Axis(fig[0, 4])
    simplify!(ax)
    colsize!(fig.layout, 4, Relative(0.04))
    # Zoom everything together
    linkaxes!(all_axs...)
    fig
end

fig

save(joinpath(basepath, "images/sim.png"), output.fig)
