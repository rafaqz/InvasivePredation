using LandscapeChange
using Unitful

include("load_data.jl")
(; cat_mass_preference, rodent_stats, params, norway_rat_studies) = 
    optimize_predation_rates_from_literature()

# Hanski D parameters taken from cat mass preference model
# These are relative, absolute value doesn't matter
predation_rates = values(map(s -> s.predation_rate, rodent_stats))
Ds = predation_rates ./ predation_rates[1]
ks = Tuple(rodent_carrycap .* u"ha")

Ns = Tuple(rodent_carrycap .* u"ha") ./ 5

# norway_rat : black_rat : mouse
# Roughly ~ 1:3:10
#
# rodents are ~80% of cat diet and they do not switch to natives much (Harper thesis etc)

cs = Tuple(max_yield_fraction) ./ t
ys = max_yield_fraction ./ t

α12 = 0.5
α13 = 0.1
α21 = 2.0
α23 = 0.2
α31 = 4.0
α32 = 3.0
# αs = (x -> 1.0).(αs)
αs = @SMatrix [0.0 α12 α13; α21 0.0 α13; α31 α32 0.0]

# cs = individuals_per_cat

hunted_rodent_mass = NamedVector(norway_rat=110, black_rat=90, mouse=15) .* u"g"
assimilated_energy_per_individual = hunted_rodent_mass .* rodent_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual
mean_prey_size = 60u"g" 
mean_prey_n = cat_energy_intake / (mean_prey_size * rodent_energy_content * assimilation_efficiency)

P = 0.00001
cs = individuals_per_cat
ys = max_yield_fraction ./ t
v = eachrow(pred_df)[1].rmax * u"yr^-1"
q = eachrow(pred_df)[1].carrycap
e = cat_energy_intake
Es = assimilated_energy_per_individual
Ncrit = 4 / u"d"
d_high = 0.2 / t
P = 0.01
Ns1 = 2 .* Ns
Ns = Tuple(rodent_carrycap .* u"ha")

for i in 1:100
    Ns = hanski_multi(P, Ns, Ds, Es, ys, αs, ks, cs, rt, d_high, t)
    P = hanski_pred(P, v, q, e, d_high, Ns, ys, Es, Ds, t)
    P
    @show (Ns, P)
end

simple_growth(N, k, r, t) = (N .* k) ./ (N .+ (k.- N) .* exp.(.-(r * t)))

N = 10
simple_growth(P, q, v, t)

