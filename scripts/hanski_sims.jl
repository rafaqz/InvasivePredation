sing Revise
using LandscapeChange
using Unitful
using StaticArrays
using StatsBase
using GLMakie
using Dispersal
using DynamicGrids
using Setfield
using Rasters
using NCDatasets
using NeutralLandscapes
using Colors
using ColorVectorSpace
GLMakie.activate!()

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
#

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
    native=8000u"kJ*d^-1",
    cleared=8000u"kJ*d^-1",
    urban=20000u"kJ*d^-1",
)

lc_exploitation_capacity = NamedVector(;
    native=NamedVector(;
        norway_rat = 0.3,
        black_rat  = 0.8,
        mouse      = 0.6,
    ),
    cleared=NamedVector(;
        norway_rat = 0.6,
        black_rat  = 0.5,
        mouse      = 0.6,
    ),
    urban=NamedVector(;
        norway_rat = 0.9,
        black_rat  = 0.6,
        mouse      = 0.6,
    ),
)

diet_overlap = NamedVector(;
    norway_rat = (0.8, 0.4),
    black_rat  = (0.8, 0.4),
    mouse      = (0.4, 0.4),
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

# Init conditions
P = 0.0001
Ns = rodent.carrycap .* u"ha" ./ 2
tspan = 1:10000
P_timeline = Vector{typeof(P)}(undef, length(tspan))
Ns_timeline = Vector{typeof(Ns)}(undef, length(tspan))
stochasticity = Param(0.05; bounds=(0.0, 0.2), label="stochasticity")
lc_names = NamedVector{keys(lc_ks)}(keys(lc_ks))
ks = lc_ks.native

supplements = map(NamedVector(native=0.0, cleared=10.0, urban=2000.0), lc_names) do v, n 
    Param(v; label="$n supplement", bounds=(0.0, 5000.0))
end
supplement = supplements.urban
model = (; P, Ns, ks, Ds, Es, ys, α, cs, rs, d_high, v, e, P_timeline, Ns_timeline, t, tspan, stochasticity, supplement, cat)

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

hanski_rule = let lc=Aux{:lc}(), lc_ks=stripparams(lc_ks), supplements=supplements, model=model
    Cell{Tuple{:cats,:rodents}}() do data, (P, Ns), I
        lc_i = get(data, lc, I)
        lc_i == 0 && return zero(P), zero(Ns)
        ks = lc_ks[lc_i]
        Ns1 = hanski_prey_timestep(P, Ns, ks, model.α, model)::typeof(Ns)
        supplement = supplements[lc_i]
        P1 = hanski_predator_timestep(P, Ns, supplement, model)::typeof(P)
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

s = 100
rodents = Raster(fill(lc_ks.native ./ 2, X(100), Y(100)))
cats = map(_ -> 0.01, rodents)
tspan = 1:1000
quantiles = (
    (0.10, 0.20),
    (0.10, 0.90),
    (0.80, 0.90),
)
lcs = map(quantiles) do (qa, qb)
    mpd = Raster(rand(MidpointDisplacement(0.35), dims(rodents)))
    a = quantile(vec(mpd), qa)
    b = quantile(vec(mpd), qb)
    classes = (<(a) => 1, a .. b => 2, >(b) => 3)
    Rasters.classify(mpd, classes; others=0, missingval=0)
end
tspan = 1:1000
init = (; cats=parent(cats), rodents=parent(rodents))
lc_mode = 2

output = ResultOutput(init; tspan=1:10, aux=(; lc=parent(lcs[lc_mode])), boundary=Wrap())
sim!(output, ruleset; printframe=true)
# Makie.heatmap(output[end].cats)
# tspan = 1800.0:t/u"yr":2000.0

output = MakieOutput(init; 
    ruleset, tspan, fps=100, 
    printframe=true, 
    aux=(; lc=parent(lcs[lc_mode])), 
    store=true, 
    ncolumns=2,
) do x
    colors = (RGBA(0.5, 0.5, 0.0), RGBA(0.5, 0.0, 0.5), RGBA(0.0, 0.5, 0.5))
    inds = [(1, 1), (2, 1), (2, 2), (2, 3)]
    I = inds[1]
    cat_ax = Axis(x.layout[1, 2]; xlabel="cat")
    c = image!(cat_ax, x.frame.cats; colormap=:magma, interpolate=false, colorrange=(0.0, 0.05))
    Colorbar(x.layout[1, 2, Right()], c; flipaxis=false)
    lc_ax = Axis(x.layout[1, 1]; xlabel="land cover")
    image!(lc_ax, lcs[lc_mode]; colormap=:viridis, interpolate=false, colorrange=(0, 4))
    rodent_axs = map(1:length(lc_ks[1]), propertynames(lc_ks[1]), maximum(lc_ks)) do i, name, k
        I = inds[i + 1]
        rodent_ax = Axis(x.layout[I...]; xlabel=string(name))
        A = rebuild(getindex.(rodents, i); missingval=0.0)
        rodent_obs = Observable{Raster}(A)
        on(x.frame.rodents) do rodents
            rodent_obs[] .= getindex.(rodents, i)
            notify(rodent_obs)
        end
        r = image!(rodent_ax, rodent_obs; colormap=:viridis, interpolate=false, colorrange=(0.0, k))
        Colorbar(x.layout[I..., Right()], r; flipaxis=false)
        rodent_ax
    end
    linkaxes!(cat_ax, lc_ax, rodent_axs...)
    nothing
end

fig = let 
    fig = Figure()
    for j in 1:3
        output = ResultOutput(init; tspan=1:500, aux=(; lc=parent(lcs[j])), boundary=Wrap())
        sim!(output, ruleset; printframe=true)
        lc_ax = Axis(fig[j, 1]; xlabel="land cover")
        image!(lc_ax, lcs[j]; colormap=:terrain, interpolate=false, colorrange=(0, 4))
        cat_ax = Axis(fig[j, 2]; xlabel="cat")
        c = image!(cat_ax, output[end].cats; colormap=:magma, interpolate=false, colorrange=(0.0, 0.05))
        Colorbar(fig[j, 2, Right()], c; flipaxis=false)
        rodent_axs = map(1:length(lc_ks[1]), propertynames(lc_ks[1]), maximum(lc_ks)) do i, name, k
            I = j, 2 + i
            rodent_ax = Axis(fig[I...]; xlabel=string(name))
            A = getindex.(output[end].rodents, i)
            r = image!(rodent_ax, A; colormap=:viridis, interpolate=false, colorrange=(0.0, k))
            Colorbar(fig[I..., Right()], r; flipaxis=false)
            rodent_ax
        end
        linkaxes!(cat_ax, lc_ax, rodent_axs...)
    end
    fig
end
fig

save(joinpath(basepath, "images/sim.png"), output.fig)
