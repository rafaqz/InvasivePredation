using Revise
using LandscapeChange
using Unitful
using StaticArrays
using GLMakie
using InvasivePredation

basepath = InvasivePredation.basepath

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
t = 1u"d"
Ds = predation_rates ./ predation_rates[1]
cs = individuals_per_cat # Base intake reequirement
# ys = map(x -> Param(x; label="yield"), max_yield_fraction) # Prey yields
ys = max_yield_fraction # Prey yields
ks = map(propertynames(rodent_carrycap), values(rodent_carrycap .* u"ha")) do k, v
    Param(v; label="$k carrycap")
end
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

param_df = DataFrame(mm)
alphapath = joinpath(basepath, "tables/alpha")
class = :native
class = :urban
class = :abandoned
class = :cleared
CSV.write(joinpath(alphapath, "$(class)_params.csv"), param_df[!, [:val, :label]])

classes = (:native, :abandoned, :cleared, :urban)
lc_α = map(classes) do class
    param_df = CSV.read(joinpath(alphapath, "$(class)_params.csv"), DataFrame)
    m = Model(α)
    m[:val] = param_df.val
    parent(m)
end


# Spatial

using DynamicGrids, Setfield

hanski_rule = let landcover = Aux{:landcover}(), lc_α=lc_α, model=model
    Cell{Tuple{:cats,:rodents}}() do data, (P, Ns), I
        lc = get(data, landcover, I)
        α = lc_α[lc]
        Ns1 = hanski_prey_timestep(P, Ns, ks, α, model1)
        P1 = hanski_predator_timestep(P, Ns, model1)
        return P1, Ns1 
    end
end
