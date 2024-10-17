using Revise
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
(; max_yield_fraction) = InvasivePredation.get_max_yield_fraction()
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

# Hanski D parameters taken from cat mass preference model
# These are relative, absolute value doesn't matter

# From other script: make this a function

t = u"yr" / 12
predation_rates = values(map(s -> s.predation_rate, rodent_stats))
Ds = predation_rates ./ predation_rates[1]
cs = individuals_per_cat # Base intake reequirement
# ys = map(x -> Param(x; label="yield"), max_yield_fraction) # Prey yields
ys = max_yield_fraction # Prey yields
# ks = map(propertynames(rodent_carrycap), values(rodent_carrycap .* u"ha")) do k, v
#     Param(v; label="$k carrycap")
# end |> NamedVector{(:norway_rat, :black_rat, :mouse)}
ks = rodent.carrycap .* u"ha"
v = cat.rmax * u"yr^-1" # Predator rmax
# q = eachrow(pred_df)[1].carrycap # Predator carrycap
e = cat.energy_intake
Es = assimilated_energy
rs = rodent.rmax
d_high = 0.2 / t
# Ncrit = 4 / u"d"
# Ns1 = 2 .* Ns

# How can we be more systematic about these
α12 = Param(1.0; label="α R.n R.r", bounds=(0.0, 2.0))
α13 = Param(1.5; label="α R.n M.m", bounds=(0.0, 2.0))
α21 = Param(0.3; label="α R.r R.n", bounds=(0.0, 2.0))
α23 = Param(1.0; label="α R.r M.m", bounds=(0.0, 2.0))
α31 = Param(0.3; label="α M.m R.n", bounds=(0.0, 2.0))
α32 = Param(0.2; label="α M.m R.r", bounds=(0.0, 2.0))
α = (
    norway_rat=(α21, α31),
    black_rat=(α12, α32),
    mouse=(α13, α23),
)

# Init conditions
P = 0.0001
Ns = rodent.carrycap .* u"ha" ./ 2
tspan = 1:10000
P_timeline = Vector{typeof(P)}(undef, length(tspan))
Ns_timeline = Vector{typeof(Ns)}(undef, length(tspan))
stochasticity = Param(0.05; bounds=(0.0, 0.2), label="stochasticity")
model = (; P, Ns, ks, Ds, Es, ys, α, cs, rs, d_high, v, e, P_timeline, Ns_timeline, t, tspan, stochasticity)

# Live interaction to find α parameters
mm = MakieModel(model) do layout, m
    ax1 = Axis(layout[1, 1])
    ax2 = Axis(layout[2, 1])
    res = lift(pred_prey_sim, m)
    P_timeline = lift(first, res)
    Rn_timeline = lift(r -> getindex.(last(r), 1), res)
    Rr_timeline = lift(r -> getindex.(last(r), 2), res)
    Mm_timeline = lift(r -> getindex.(last(r), 3), res)
    Makie.lines!(ax1, Rn_timeline; label=rodent.labels[1])
    Makie.lines!(ax1, Rr_timeline; label=rodent.labels[2])
    Makie.lines!(ax1, Mm_timeline; label=rodent.labels[3])
    Makie.lines!(ax2, P_timeline; label="Cat")
    linkxaxes!(ax1, ax2)
    axislegend(ax1)
    axislegend(ax2)
end

lc_ks = (
    native    = NamedVector(norway_rat=20.0, black_rat=50.0, mouse=40.0), 
    cleared   = NamedVector(norway_rat=30.0, black_rat=10.0, mouse=50.0), 
    urban     = NamedVector(norway_rat=200.0, black_rat=150.0, mouse=200.0),
)

supplements = map(x -> x * u"kJ*d^-1", (native=0.0, cleared=0.0, urban=1.0))

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
    # maskbehavior=Dispersal.CheckMaskEdges()
)

rodent_λ = NamedVector{propertynames(Ns)}((
    Param(0.5, label="R.n λ", bounds=(0.00001, 1.0)),
    Param(0.3, label="R.r λ", bounds=(0.00001, 1.0)),
    Param(0.2, label="M.n λ", bounds=(0.00001, 1.0)),
))
ExponentialKernel([1.0, 2.0])(10.0)

rodent_spread_rule = OutwardsDispersal{:rodents}(; 
    stencil=Moore(4),
    formulation=ExponentialKernel(rodent_λ),
    maskbehavior=Dispersal.CheckMaskEdges()
)

rodent_allee_rule = AlleeExtinction{:rodents}(NamedVector(; norway_rat=0.1, black_rat=0.4, mouse=0.5))

ruleset = Ruleset(hanski_rule, rodent_spread_rule, cat_spread_rule, rodent_allee_rule)

# rodents = map(_ -> ks ./ 2, lc[Ti=1])
rodents = Raster(fill(ks ./ 2, X(64), Y(64)))
cats = map(_ -> 0.01, rodents)
tspan = 1:1000
lcs = map(0:0.1:0.2) do x
    a = 0.35 + x
    b = 0.55 + x
    classes = >(b) => 1, a .. b => 2, <=(a) => 3
    mpd = Raster(rand(MidpointDisplacement(0.4), dims(rodents)))
    Rasters.classify(mpd, classes; others=0, missingval=0)
end
tspan=1:1000
lc_mode = 1
init = (; cats=parent(cats), rodents=parent(rodents))

output = ResultOutput(init; tspan, aux=(; lc=parent(lcs[lc_mode])), boundary=Wrap())
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
    Colorbar(x.layout[I..., Right()], c; flipaxis=false)
    lc_ax = Axis(x.layout[1, 1]; xlabel="land cover")
    image!(lc_ax, lcs[lc_mode]; colormap=:viridis, interpolate=false, colorrange=(0, 4))
    # rodent_axs = map(1:length(ks), propertynames(ks), ks) do i, name, k
    #     I = inds[i + 1]
    #     rodent_ax = Axis(x.layout[I...]; xlabel=string(name))
    #     A = rebuild(getindex.(rodents, i); missingval=0.0)
    #     rodent_obs = Observable{Raster}(A)
    #     on(x.frame.rodents) do rodents
    #         rodent_obs[] .= getindex.(rodents, i)
    #         notify(rodent_obs)
    #     end
    #     r = image!(rodent_ax, rodent_obs; colormap=:viridis, interpolate=false, colorrange=(0.0, 2k))
    #     Colorbar(x.layout[I..., Right()], r; flipaxis=false)
    #     rodent_ax
    # end
    maxs = Ref(ks)
    rodent_ax = Axis(x.layout[1, 3]; xlabel="Rodents")
    rgb(x) = sum(min.(1.0, x ./ maxs[]) .* colors)
    A = rebuild(map(rgb, rodents); missingval=RGBA(0.0, 0.0, 0.0, 0.0))
    rodent_obs = Observable{Raster}(A)
    on(x.frame.rodents) do rodents
        maxs[] = max.(maximum(rodents), maxs[])
        rodent_obs[] .= rgb.(rodents)
        notify(rodent_obs)
    end
    r = image!(rodent_ax, rodent_obs; interpolate=false)
    # Colorbar(x.layout[I..., Right()], r; flipaxis=false)
    linkaxes!(cat_ax, lc_ax, rodent_ax)
    nothing
end

# save(joinpath(basepath, "images/sim.png"), output.fig)
