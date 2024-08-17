

#=
Another paper using this simulation:

Mesopredator release of rats can only happen when there is a high-density
primary food-source for cats *other than rodents that compete with each other*.
This is likely limited to urban feeding/refuse,
high density sea bird (and similar) colonies, or rabbits.

Demonstration: island with or without rabbits and only rats
- grasslands with rabbits will show meopredator release
- forests will not show mesopredator release
=#

# using DynamicGrids, Dispersal

# space = 8, 8
# bird_colony = rand(space...) .> 0.95
# landcover = rand(1:3, space...)
# rodents = fill(rodent_carrycap .* u"ha", space...)
# cats = fill(0.01, space...)
# init = (; cats, rodents)
# tspan = 1:100
# rodent_growth = LogisticGrowth(; carrycap=rodent_carrycap * u"ha", timestep=1u"yr")
# cat_growth = let max_yield_fraction=max_yield_fraction
#     Cell{Tuple{:cat,:rodents}}() do data, (cat, rodents), I
#         max_yield = rodent .* max_yield_fraction
#         max_supported_cats = sum(net_energy_per_rodent .*  max_yield)
#         if max_supported_cats < cat
#             (N * k) / (N + (k - N) * exp(-rt))
#         else
#             N * exp(rt)
#         end
#         reduced_rodents = rodents .- take
#         (grown_cat, reduced_rodents)
#     end
# end

# ruleset = Ruleset(cat_growth, rodent_growth; boundary=Wrap())

# output = MakieOutput(init;
#     ruleset,
#     tspan,
#     aux = (; black_rat_predation_rate, landcover),
# )

# R = NamedVector{Tuple(rodent_names),3}
# Hanski D parameters
# These are relative, absolute value doesn't matter
Ds = predation_rates ./ predation_rates[1]

# norway_rat : black_rat : mouse
# Roughly ~ 1:3:10
#
# rodents are ~80% of cat diet and they do not switch to natives much (Harper thesis etc)

rodent_keys = Tuple(rodent_names)
# pred_funcs = (;
#     black_rat  = p -> -0.1f0p.norway_rat - 0.1f0p.mouse + 0.5f0p.native + 0.3f0p.abandoned + 0.3f0p.forestry + 1p.urban,
#     norway_rat = p -> -0.1f0p.black_rat - 0.1f0p.mouse + 1.5f0p.urban - 0.2f0p.native,
#     mouse =      p -> -0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
#     pig =        p -> 0.0f0p.native - 0.0f3p.abandoned - 2f0p.urban - 1.0f0p.cleared,
#     wolf_snake = p -> -0.2f0p.black_rat + 0.3f0p.mouse - 0.5f0p.urban + 0.3f0p.native,
#     macaque =    p -> 1.0f0p.abandoned + 0.7f0p.forestry + 0.4f0p.native - 1.0f0p.urban - 0.8f0p.cleared
# )[rodent_keys]

# carrycap = NamedVector(;
#     cat =        0.02,
#     black_rat =  30.0, # This may be up to 100/ha? See Harper & Bunbury 2015.
#     norway_rat = 15.0, # This one is more of a guess
#     mouse =      52.0,
#     pig =        0.02, # Tuned to have ~600 total in mauritius in 2009 (actually 714, more significant digits are not justifiable).
#     wolf_snake = 10.0, # As is this one
#     macaque =    0.5,
# )[rodent_keys]

# populations = carrycap ./ 2
# local_inputs = NamedVector(native=1.0, cleared=0.0, abandoned=0.0, urban=0.0, forestry=0.0)
# calc_carrycaps(local_inputs, populations, carrycap, pred_funcs)

