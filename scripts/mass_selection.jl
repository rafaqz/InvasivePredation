using Revise
using LandscapeChange
using Unitful
using InvasivePredation
using StaticArrays
using GLMakie

include("load_data.jl")
(; cat_mass_preference, rodent_stats, params, norway_rat_studies) = 
    optimize_predation_rates_from_literature()

# Hanski D parameters taken from cat mass preference model
# These are relative, absolute value doesn't matter

# From other script: make this a function
max_yield_fraction = (0.10578653307702114, 0.1344295847234702, 0.17522617562570525) .* u"yr^-1"
predation_rates = values(map(s -> s.predation_rate, rodent_stats))

# How can we be more systematic about these
α12 = 0.1
α13 = 0.1
α21 = 1.0
α23 = 0.1
α31 = 1.3
α32 = 1.3
αs = @SMatrix [0.0 α12 α13; α21 0.0 α13; α31 α32 0.0]
# αs = (x -> 1.0).(αs)

# cs = individuals_per_cat

hunted_rodent_mass = NamedVector(norway_rat=110, black_rat=90, mouse=15) .* u"g"
assimilated_energy_per_individual = hunted_rodent_mass .* rodent_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual
mean_prey_size = 60u"g" 
mean_prey_n = cat_energy_intake / (mean_prey_size * rodent_energy_content * assimilation_efficiency)

t = 1u"d"
Ds = predation_rates ./ predation_rates[1]
cs = individuals_per_cat # Base intake reequirement
ys = max_yield_fraction # Prey yields
v = eachrow(pred_df)[1].rmax * u"yr^-1" # Predator rmax
q = eachrow(pred_df)[1].carrycap # Predator carrycap
e = cat_energy_intake
Es = assimilated_energy_per_individual
rs = rodent_rmax 
# Ncrit = 4 / u"d"
d_high = 0.2 / t
# Ns1 = 2 .* Ns

P = 0.0001 
Ns = Tuple(rodent_carrycap .* u"ha")
P_timeline = [P]
Ns_timeline = [Ns]
for i in 1:10000
    ks1 = if i > 3000 && i < 3100
        100 .* ks # mast year
    else
        ks
    end
    Ns = hanski_multi(P, Ns, Ds, Es, ys, αs, ks1, cs, rs, d_high, t)
    P = hanski_pred(P, v, e, d_high, Ns, ys, Es, Ds, t)
    push!(P_timeline, P)
    push!(Ns_timeline, Ns)
    @show (Ns, P)
end


fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
Makie.lines!(ax1, map(first, Ns_timeline); label=rodent_labels[1])
Makie.lines!(ax1, map(x -> x[2], Ns_timeline); label=rodent_labels[2])
Makie.lines!(ax1, map(last, Ns_timeline); label=rodent_labels[3])
Makie.lines!(ax2, P_timeline; label="cat")
axislegend(ax1)


# simple_growth(N, k, r, t) = (N .* k) ./ (N .+ (k.- N) .* exp.(.-(r * t)))
# N = 10
# simple_growth(P, q, v, t)

