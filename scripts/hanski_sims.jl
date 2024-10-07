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
GLMakie.activate!()

using InvasivePredation

basepath = InvasivePredation.basepath
alphapath = joinpath(basepath, "tables/alpha")
(; cat, rodent) = s = InvasivePredation.load_settings()
(; cat_mass_preference, rodent_stats, norway_rat_params, norway_rat_studies) =
    fit_distributions_to_literature()
(; max_yield_fraction) = InvasivePredation.get_max_yield_fraction()
(; assimilated_energy, individuals_per_cat) = InvasivePredation.get_cat_energetics(cat, rodent, rodent_stats)

lc_path = joinpath(basepath, "data/lc_predictions_mus.nc")
lc = RasterStack(lc_path) |>
    x -> Rasters.maybeshiftlocus(Rasters.Start(), x) |>
    x -> Rasters.set(x, Ti => Int.(lookup(x, Ti))) |>
    x -> Rasters.aggregate(Rasters.Start(), x, (X(aggfactor), Y(aggfactor))) |>
    x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |> NamedTuple
lc = lc.native .* 1 .+ lc.cleared .* 2 .+ lc.abandoned .* 3  .+ lc.forestry .* 3 .+ lc.urban .* 4
m = lc .!= 0

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
α12 = Param(0.6; label="α R.n R.r", bounds=(0.0, 2.0))
α13 = Param(1.0; label="α R.n M.m", bounds=(0.0, 2.0))
α21 = Param(0.4; label="α R.r R.n", bounds=(0.0, 2.0))
α23 = Param(0.8; label="α R.r M.m", bounds=(0.0, 2.0))
α31 = Param(0.4; label="α M.m R.n", bounds=(0.0, 2.0))
α32 = Param(0.5; label="α M.m R.r", bounds=(0.0, 2.0))
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
model = (; P, Ns, ks, Ds, Es, ys, α, cs, rs, d_high, v, e, P_timeline, Ns_timeline, t, tspan)

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

# Manual parameter definition

# class = :native
# class = :urban
# class = :abandoned
# class = :cleared

# CSV.write(joinpath(alphapath, "$(class)_params.csv"), DataFrame(mm)[!, [:val, :label]])

classes = (:native, :abandoned, :cleared, :urban)
# lc_α = map(classes) do class
#     param_df = CSV.read(joinpath(alphapath, "$(class)_params.csv"), DataFrame)
#     m = Model(α)
#     m[:val] = param_df.val
#     parent(m)
# end

# Spatial

# ks = rodent_carrycap .* u"ha"

aggfactor = 4

lc_ks = (
    native    = NamedVector(norway_rat=80.0, black_rat=90.0, mouse=150.0), 
    cleared   = NamedVector(norway_rat=30.0, black_rat= 30.0, mouse=100.0), 
    abandoned = NamedVector(norway_rat=30.0, black_rat= 30.0, mouse=100.0), 
    urban     = NamedVector(norway_rat=100.0, black_rat=200.0, mouse=350.0),
)

hanski_rule = let lc=Aux{:lc}(), lc_ks=stripparams(lc_ks), model=model
    Cell{Tuple{:cats,:rodents}}() do data, (P, Ns), I
        lc_i = get(data, lc, I)
        lc_i == 0 && return zero(P), zero(Ns)
        # α = lc_α[lc_i]
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
    Param(0.5, label="R.n λ", bounds=(0.00001, 4.0)),
    Param(0.3, label="R.r λ", bounds=(0.00001, 4.0)),
    Param(0.2, label="M.n λ", bounds=(0.00001, 4.0)),
))
ExponentialKernel([1.0, 2.0])(10.0)

rodent_spread_rule = OutwardsDispersal{:rodents}(; 
    stencil=Moore(4),
    formulation=ExponentialKernel(rodent_λ),
    maskbehavior=Dispersal.CheckMaskEdges()
)

rodent_allee_rule = AlleeExtinction{:rodents}(1.0)

ruleset = Ruleset(hanski_rule, rodent_spread_rule, cat_spread_rule, rodent_allee_rule)

rodents = map(lc[Ti=1]) do _
    ks ./ 2
end
# rodents[X=1:100] .= (zero(ks),)
cats = map(lc[Ti=1]) do _
    0.01
end
init = (; cats, rodents)

# ruleset = Ruleset(hanski_rule)
tspan = 1800.0:t/u"yr":2000.0

output = ResultOutput(init; tspan, aux=(; lc))
# sim!(output, ruleset; printframe=true)
# @profview sim!(output, ruleset; printframe=true)
init = (; a=rand(100, 100))
tspan=1:10
# TODO use observables for this...
i = Ref(0)
output = TransformedOutput(init; tspan, eltype=Union{Missing,NamedTuple}) do data
    i[] += 1
    i[] in (2, 10, 12) ? data : missing
end

output = MakieOutput(init; ruleset, tspan, fps=100, printframe=true, aux=(; lc), mask=m[Ti=1], store=true) do x
    inds = [(1, 1), (2, 1), (1, 2), (2, 2)]
    I = inds[1]
    cat_ax = Axis(x.layout[I...]; xlabel="cat")
    c = image!(cat_ax, x.frame.cats; colormap=:magma, interpolate=false, colorrange=(0.0, 0.05))
    Colorbar(x.layout[I..., Right()], c; flipaxis=false)
    rodent_axs = map(1:length(ks), propertynames(ks), ks) do i, name, k
        I = inds[i + 1]
        rodent_ax = Axis(x.layout[I...]; xlabel=string(name))
        A = rebuild(getindex.(rodents, i); missingval=0.0)
        rodent_obs = Observable{Raster}(A)
        on(x.frame.rodents) do rodents
            rodent_obs[] .= getindex.(rodents, i)
            notify(rodent_obs)
        end
        r = image!(rodent_ax, rodent_obs; colormap=:viridis, interpolate=false, colorrange=(0.0, 2k))
        Colorbar(x.layout[I..., Right()], r; flipaxis=false)
        rodent_ax
    end
    linkaxes!(cat_ax, rodent_axs...)
    nothing
end

# save(joinpath(basepath, "images/sim.png"), output.fig)
