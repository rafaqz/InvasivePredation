using Revise
using LandscapeChange
using Unitful
using StaticArrays
using StatsBase
using GLMakie
using InvasivePredation
using Dispersal

basepath = InvasivePredation.basepath
alphapath = joinpath(basepath, "tables/alpha")

function pred_prey_sim(model)
    (; P, Ns, P_timeline, Ns_timeline, ks, α, tspan) = model
    P1 = P
    Ns1 = Ns
    @inbounds for i in tspan
        # ks1 = if i > 3200 && i < 3200 || i > 6200 && i < 6200
        #     2 .* ks # mast year
        # else
        #     ks
        # end
        Ns2 = hanski_prey_timestep(P1, Ns1, ks, α, model)::typeof(Ns)
        P1 = hanski_predator_timestep(P1, Ns1, model)::typeof(P)
        Ns1 = Ns2
        P_timeline[i] = P1
        Ns_timeline[i] = Ns1
    end
    return P_timeline, Ns_timeline
end

include("load_settings.jl")
(; cat_mass_preference, rodent_stats, norway_rat_params, norway_rat_studies) =
    fit_distributions_to_literature()

# Hanski D parameters taken from cat mass preference model
# These are relative, absolute value doesn't matter

# From other script: make this a function
max_yield_fraction = (0.10578653307702114, 0.1344295847234702, 0.17522617562570525) .* u"yr^-1" .* 12
predation_rates = values(map(s -> s.predation_rate, rodent_stats))
hunted_rodent_mass = NamedVector(norway_rat=110, black_rat=90, mouse=15) .* u"g"
assimilated_energy_per_individual = hunted_rodent_mass .* rodent_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual
mean_prey_size = 60u"g"
mean_prey_n = cat_energy_intake / (mean_prey_size * rodent_energy_content * assimilation_efficiency)
t = u"yr" / 12
Ds = predation_rates ./ predation_rates[1]
cs = individuals_per_cat # Base intake reequirement
# ys = map(x -> Param(x; label="yield"), max_yield_fraction) # Prey yields
ys = max_yield_fraction # Prey yields
ks = map(propertynames(rodent_carrycap), values(rodent_carrycap .* u"ha")) do k, v
    Param(v; label="$k carrycap")
end
ks = rodent_carrycap .* u"ha"
v = eachrow(pred_df)[1].rmax * u"yr^-1" # Predator rmax
# q = eachrow(pred_df)[1].carrycap # Predator carrycap
e = cat_energy_intake
Es = assimilated_energy_per_individual
rs = rodent_rmax
d_high = 0.2 / t
# Ncrit = 4 / u"d"
# Ns1 = 2 .* Ns

# How can we be more systematic about these
α12 = Param(0.7; label="α R.n R.r")
α13 = Param(0.3; label="α R.n M.m")
α21 = Param(1.3; label="α R.r R.n")
α23 = Param(1.0; label="α R.r M.m")
α31 = Param(1.4; label="α M.m R.n")
α32 = Param(1.0; label="α M.m R.r")
α = @SMatrix [0.0 α12 α13
              α21 0.0 α23
              α31 α32 0.0]

# Init conditions
P = 0.0001
Ns = rodent_carrycap .* u"ha" ./ 2
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
    Makie.lines!(ax1, Rn_timeline; label=rodent_labels[1])
    Makie.lines!(ax1, Rr_timeline; label=rodent_labels[2])
    Makie.lines!(ax1, Mm_timeline; label=rodent_labels[3])
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
lc_α = map(classes) do class
    param_df = CSV.read(joinpath(alphapath, "$(class)_params.csv"), DataFrame)
    m = Model(α)
    m[:val] = param_df.val
    parent(m)
end

# Spatial

using DynamicGrids
using Setfield
using Rasters
using NCDatasets

aggfactor = 4

hanski_rule = let lc = Aux{:lc}(), lc_α=stripparams(lc_α), model=model, ks=ks
    Cell{Tuple{:cats,:rodents}}() do data, (P, Ns), I
        lc_i = get(data, lc, I)
        lc_i == 0 && return zero(P), zero(Ns)
        α = lc_α[lc_i]
        Ns1 = hanski_prey_timestep(P, Ns, ks, α, model)
        P1 = hanski_predator_timestep(P, Ns, model)
        return P1, Ns1
    end
end

lc_path = joinpath(basepath, "data/lc_predictions_mus.nc")
lc = RasterStack(lc_path) |>
    x -> Rasters.maybeshiftlocus(Rasters.Start(), x) |>
    x -> Rasters.set(x, Ti => Int.(lookup(x, Ti))) |>
    x -> Rasters.aggregate(Rasters.Start(), x, (X(aggfactor), Y(aggfactor))) |>
    x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |> NamedTuple
lc = lc.native .* 1 .+ lc.cleared .* 2 .+ lc.abandoned .* 3  .+ lc.forestry .* 3 .+ lc.urban .* 4
ks = rodent_carrycap .* u"ha"

rodents = map(lc[Ti=1]) do _
    ks / 2
end
cats = map(lc[Ti=1]) do _
    0.01
end
init = (; cats, rodents)

kern = DispersalKernel(
    stencil=Moore(2),
    formulation=ExponentialKernel(Param(1.0, label="λ", bounds=(0.0, 2.0))),
)

pred_spread_rule = let aggfactor=aggfactor, kern=kern
    SetNeighbors{:cats}(kern) do data, hood, N, I
        N === zero(N) && return nothing
        cellsize = (100 * 4)
        sum = zero(N)
        # # Randomise hood starting position to avoid directional artifacts in output
        start = rand(0:length(hood)-1)
        @inbounds for ix in eachindex(hood)
            # Rotate indices in relation to starting point
            i = start + ix
            if i > length(hood)
                i = i - length(hood)
            end
            Ih = DynamicGrids.indices(hood, I)[i]
            # d = DynamicGrids.distances(hood)[i]
            propagules = trunc(N * rand(typeof(N)) ^ 3)
            sum1 = sum .+ propagules
            if sum1 > N
                propagules = min(sum1, n) - sum
                sum = sum + propagules
            else
                sum = sum1
            end
            add!(data[:cats], propagules, Ih...)
        end
        @inbounds sub!(data[:cats], sum, I...)
        return nothing
    end
end # let

rodent_spread_rule = InwardsDispersal{:rodents}(; 
    stencil=Moore(1),
    formulation=ExponentialKernel(Param(0.1, label="λ", bounds=(0.0, 2.0))),
)

ruleset = Ruleset(hanski_rule, pred_spread_rule)#, rodent_spread_rule)
tspan = 1900.0:t/u"yr":2000+1

output = ResultOutput(init; tspan, aux=(; lc))
sim!(output, ruleset; printframe=true)

output = MakieOutput(init; ruleset, tspan, printframe=true, aux=(; lc)) do x
    inds = [(1, 1), (2, 1), (1, 2), (2, 2)]
    cat_ax = Axis(x.layout[inds[1]...]; xlabel="cat")
    image!(cat_ax, x.frame.cats; colormap=:magma, interpolate=false, colorrange=(0.0, 0.03))
    rodent_axs = map(1:length(ks), propertynames(ks), ks) do i, name, k
        rodent_ax = Axis(x.layout[inds[i + 1]...]; xlabel=string(name))
        rodent_obs = Observable{Raster}(rebuild(getindex.(rodents, i); missingval=0.0))
        on(x.frame.rodents) do rodents
            rodent_obs[] .= getindex.(rodents, i)
            notify(rodent_obs)
        end
        image!(rodent_ax, rodent_obs; colormap=:viridis, interpolate=false, colorrange=(0.0, k))
        rodent_ax
    end
    linkaxes!(cat_ax, rodent_axs...)
    nothing
end

save(joinpath(basepath, "images/sim.png"), output.fig)
