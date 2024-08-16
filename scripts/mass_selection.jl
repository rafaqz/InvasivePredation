deleteat!(Base.LOAD_PATH, 2:3)

using CSV
using DataFrames
using Distributions
using LandscapeChange
using Optimization
using OptimizationOptimJL
using StaticArrays
using Unitful
using GLMakie
# using CairoMakie
# CairoMakie.activate!()

seasonal_k(k, seasonality, nsteps, s) = k + k * seasonality * sin(2π * s / nsteps)

function rodent_func(x, p)
    (; k, rt, nsteps, years, seasonality, fixed_take, replicates) = p
    replicates = map(1:replicates) do i
        N = k
        # We always have to leave some population or it cant recover.
        minleft = oneunit(N) * 0.01
        total_taken = zero(N)

        # Burn in
        for _ in 1:5
            for s in 1:nsteps
                # @show N
                k_season = seasonal_k(k, seasonality, nsteps, s)
                yield = rand(LogNormal(log(x[1]), p.std))
                N = (N .* k_season) ./ (N .+ (k_season .- N) .* exp.(.-(rt / (nsteps * u"yr^-1"))))
                # println()
                # @show k_season N
                # Cant take the whole population
                taken = if fixed_take
                    max(-N + minleft, k * -yield)
                else
                    max(-N + minleft, N * -yield)
                end
                N += taken
            end
        end
        total_taken = zero(N)
        for _ in 1:years
            for s in 1:nsteps
                # @show N
                k_season = seasonal_k(k, seasonality, nsteps, s)
                yield = rand(LogNormal(log(x[1]), p.std))
                N = (N .* k_season) ./ (N .+ (k_season .- N) .* exp.(.-(rt / (nsteps * u"yr^-1"))))
                # println()
                # @show k_season N
                # Cant take the whole population
                taken = if fixed_take
                    max(-N + minleft, k * -yield)
                else
                    max(-N + minleft, N * -yield)
                end
                @assert taken <= zero(taken) "taken: $taken N: $N yield: $yield"
                @assert taken >= -N "taken: $taken N: $N yield: $yield"
                N += taken
                total_taken += taken
            end
        end
        return (taken=total_taken / years, final_pop=N)
    end
    return (taken=mean(first, replicates), final_pop=mean(last, replicates))
end
rodent_take_unitless(x, p) = ustrip(rodent_take(x, p))
rodent_take(x, p) = rodent_func(x, p)[1]
rodent_pop(x, p) = rodent_func(x, p)[2]

pred_names = ["cat", "black_rat", "norway_rat", "mouse", "pig"]
# run(`libreoffice tables/animals.csv`)
pred_df = filter(CSV.read("tables/animals.csv", DataFrame)) do row
    row.name in pred_names
end

cat_data = map(identity, NamedTuple(filter(r -> r.name == "cat", pred_df)[1, :]))
rodent_df = pred_df[2:4, :]
rodent_names = map(Symbol, rodent_df.name)
rodent_labels = titlecase.(replace.(rodent_df.name, ("_" => " ",)))
R = NamedVector{Tuple(rodent_names),3}

cat_energy_intake = 2131u"kJ/d"
assimilation_efficiency = 0.84
rodent_energy_content = 6.24u"kJ/g"
# rodent_mass = R(rodent_df.mass) .* u"g"
# Use mean from actual sizes taken for Norway rats, and guess the others
# For mouse we just use a size 15, approximately half mean (men and mice paper)
rodent_mass = NamedVector(norway_rat=110, black_rat=90, mouse=15) .* u"g"
assimilated_energy_per_individual = rodent_mass .* rodent_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual

rt = R(rodent_df.rmax) .* u"yr^-1"
rodent_carrycap = k = R(rodent_df.carrycap) .* u"ha^-1"

function optimise_hunting(rodent_params)
    p = first(rodent_params)
    optimal_takes = map(rodent_params) do p
        x0 = [0.1]
        prob = OptimizationProblem(rodent_take_unitless, x0, p; lb=[0.0], ub=[1.0])
        sol = solve(prob, SAMIN(); maxiters=10000)
        sol.u
    end
    optimal_pops = map(rodent_pop, optimal_takes, rodent_params)
    optimal_caught = map(rodent_take, optimal_takes, rodent_params) ./ u"yr" .* -1
    max_supported_cats = uconvert.(u"km^-2", optimal_caught ./ individuals_per_cat .* p.fraction_eaten)
    # @show optimal_pops optimal_caught max_supported_cats
    takes = NamedTuple{map(k -> Symbol(:taken_, k), propertynames(optimal_takes))}(values(optimal_takes))
    cats = NamedTuple{map(k -> Symbol(:cat_, k), propertynames(max_supported_cats))}(values(max_supported_cats))
    (; std=p.std, seasonality=p.seasonality, map(first, takes)..., cats..., NamedTuple(optimal_pops)...)
end

# nsteps = 365
# 12 steps (months) is a close enough approximation of continuous
nsteps = 12
years = 10
fraction_eaten = 0.72 # McGregor 2015
seasonality = 0.5
replicates = 1

# seasonality_rates = map(0.0:0.1:0.9) do seasonality
#     rodent_params = map(R(rodent_names)) do rodent
#         (; k=k[rodent], rt=rt[rodent], fixed_take=true, std=0.1, replicates=100, seasonality, nsteps, years, fraction_eaten)
#     end
#     optimise_hunting(rodent_params)
# end |> DataFrame

# rodent_params = map(R(rodent_names)) do rodent
#     (; k=k[rodent], rt=rt[rodent], fixed_take=true, std=0.2, nsteps, years, seasonality, fraction_eaten, replicates)
# end
# p = rodent_params[3]
# x0 = 0.075
# rodent_func(x0, p)
# optimise_hunting(rodent_params)

stochastic_rates = map((0.0:0.2:1.5).^2) do std
    rodent_params = map(R(rodent_names)) do rodent
        (; k=k[rodent], rt=rt[rodent], replicates=25, fixed_take=false, nsteps, years, std, seasonality, fraction_eaten)
    end
    optimise_hunting(rodent_params)
end |> DataFrame
stochastic_rates
max_yield_fraction = Tuple(DataFrame(stochastic_rates)[1, 3:5])
# Why is this relationship slightly under 1
# max_yield_fraction ./ (takes[1] / ustrip(rt[1]))

# k does not matter - we get the same yield rate regardless of k
# so the model is general to any k


# Stochasticity reduces yield
colors = [:red, :lightblue, :yellow]
x = stochastic_rates.std
fig = Figure(; title="Effects of monthly stochasticity on optimal hunting yields", size=(600, 800));
labels = ["Mean fraction hunted per month", "Supported cats per km^2", "Rodents per hectare"]
axs = map(1:3, labels) do i, ylabel
    ax = if i == 3
        Axis(fig[i, 1];
            xlabel="Scale of log-normal hunting stochasticity", ylabel,
        )
    else
        Axis(fig[i, 1]; ylabel, xticksvisible=false, xticklabelsvisible=false)
    end
    # hidespines!(ax)
    ax
end
rodent_plots = map(1:3) do i
    color = colors[i]
    label = rodent_labels[i]
    p1 = Makie.lines!(axs[1], x, ustrip.(stochastic_rates[:, i+1]);
        label, color
    )
    p2 = Makie.lines!(axs[2], x, ustrip.(stochastic_rates[:, i+4]);
        label, color
    )
    p3 = Makie.lines!(axs[3], x, ustrip.(stochastic_rates[:, i+7]);
        label, color
    )
    (p1, p2, p3)
end
axislegend(axs[1])
display(fig)
save("images/cat_hunting_stochasticty.png", fig)


using Distributions
using NLopt, OptimizationNLopt

# Childs 1986
max_rat_size = 100:100:700
center_rat_size = 50:100:650
max_rat_25_size = 25:25:200
center_rat_25_size = 12.5:25:187.5
killed_rats_25_childs = [0, 8, 8, 3, 1, 1, 0, 1]
killed_rats_childs = [19, 4, 0, 0, 0, 0, 0]
trapped_rats_childs = [5, 5, 9, 29, 44, 14, 3]
sum(center_rat_25_size .* killed_rats_25_childs) ./ sum(killed_rats_25_childs)

# Glass 2009
trapped_rats_glass = [18, 38, 49, 167, 175, 98, 0]
killed_rats_glass = [15, 13, 4, 2, 1, 0, 0]
center_rat_size .* trapped_rats_glass
center_rat_size .* killed_rats_glass
cat_preference_glass = LogNormal(log(60), 1.1);
rat_pdfs_glass = pdf.((cat_preference_glass,), center_rat_size) .* trapped_rats_glass
killed_rats_glass ./ sum(killed_rats_glass)
rat_pdfs_glass ./ sum(rat_pdfs_glass)

mean_glass = sum(center_rat_size .* killed_rats_glass) ./ sum(killed_rats_glass)
mean_childs = sum(center_rat_25_size .* killed_rats_25_childs) ./ sum(killed_rats_25_childs)
n_glass = sum(killed_rats_glass)
n_childs = sum(killed_rats_25_childs)
n_total = n_glass + n_childs
frac_glass = n_glass / n_total
frac_childs = n_childs / n_total
mean_caught_norway_rat_mass = mean_glass .* frac_glass + mean_childs .* frac_childs

trapped_rats = trapped_rats_glass .+ trapped_rats_childs
killed_rats = killed_rats_glass .+ killed_rats_childs

function get_pref(x, p)
    μ, σ = x
    cat_preference = LogNormal(log(μ), σ);
    rat_pdfs = pdf.((cat_preference,), p.center_rat_size) .* p.trapped_rats
    killed_means = p.killed_rats ./ sum(p.killed_rats)
    predicted_means = rat_pdfs ./ sum(rat_pdfs)
    sum((killed_means .- predicted_means) .^ 2)
end

x = [90.0, 0.6]
p = (; trapped_rats, killed_rats, center_rat_size)
prob = OptimizationProblem(get_pref, x, p; lb=[0.0, 0.0], ub=[500.0, 10.0])
get_pref(x, p)
sol = solve(prob, NLopt.LN_NELDERMEAD())
x = sol.u
get_pref(x, p)
μ, σ = x
cat_preference = LogNormal(log(μ), σ)
Makie.plot(cat_preference)


# From the diet review paper
# these are means prey size accross whole studies
# using Statistics

# x0 = [20.0, -10.0]
# p = (; cat_mean_prey_sizes)
# get_pref(x0, p)
# prob = OptimizationProblem(prey_size, x0, p)
# sol = solve(prob, NelderMead())

# x0 = [0.1]
# p = (; cat_mean_prey_sizes)
# prey_size(x0, p)
# prob = OptimizationProblem(prey_size, x0, p)
# sol = solve(prob, NelderMead(); )
# Makie.hist!(cat_mean_prey_sizes)
# p = Makie.plot(cat_preference)
# xlims!(p.axis, (0, 300))

cat_mean_prey_sizes = [11.0, 33.0, 17.9, 4.3, 15.5, 23.6, 24.0, 31.2, 26.5, 23.3, 9.3, 26.2, 29.0, 21.2, 349.7, 3.64, 21.8, 35.9, 13.6, 8.2, 228.6, 249.0, 42.7, 16.5, 15.5, 60.7, 32.8, 184.9, 34.1, 36.8, 7.4, 102.0, 241.2, 38.5, 13.3, 26.4, 7.6, 24.3, 123.2, 39.2, 1.77, 76.7, 13.5, 18.3, 30.3, 34.4]
sort!(cat_mean_prey_sizes)
mean(cat_mean_prey_sizes)
cat_preference = LogNormal(log(50), 0.8);

# Wilson et al 2007
black_rat_max_mass = 40:40:240
black_rat_center_mass = 20:40:220
black_rat_trap_rate = [0.08, 0.20, 0.26, 0.31, 0.14, 0.01]
mean_black_rat_mass = sum(black_rat_trap_rate .* black_rat_center_mass)
black_rat_mass_distribution = black_rat_trap_rate .* black_rat_center_mass ./ mean_black_rat_mass

# Glass 1988
norway_rat_max_mass = 100:100:500
norway_rat_center_mass = norway_rat_max_mass .- 50
male_norway_rat_trap_rate = [0.07, 0.3, 0.43, 0.16, 0.04]
female_norway_rat_trap_rate = [0.14, 0.33, 0.34, 0.19, 0.0]
norway_rat_trap_rate = (male_norway_rat_trap_rate .+ female_norway_rat_trap_rate) ./ 2
mean_norway_rat_mass = sum(norway_rat_trap_rate .* norway_rat_center_mass)
norway_rat_mass_distribution = norway_rat_trap_rate .* norway_rat_center_mass ./ mean_norway_rat_mass

# made up from shifted rat samples
mouse_max_mass = 10:10:60
mouse_center_mass = 5:10:55
mouse_trap_rate = [0.08, 0.20, 0.26, 0.31, 0.14, 0.01]
mean_mouse_mass = sum(mouse_trap_rate .* mouse_center_mass)
mouse_mass_distribution = (mouse_trap_rate .* mouse_center_mass) ./ mean_mouse_mass
mouse_mass_yields = mouse_trap_rate .* mouse_mass_distribution

norway_rat_pdfs = pdf.((cat_preference,), norway_rat_center_mass)
norway_rat_trap_rate
norway_rat_catch_rates = norway_rat_pdfs .* norway_rat_trap_rate
norway_rat_normalised_rates = norway_rat_catch_rates ./ sum(norway_rat_trap_rate)
norway_rat_normalised_mass_yield = sum(norway_rat_normalised_rates .* norway_rat_mass_distribution .* mean_norway_rat_mass)
norway_rat_mass_yields = sum(norway_rat_catch_rates .* norway_rat_mass_distribution)
norway_rat_predation_rate = sum(norway_rat_catch_rates)

black_rat_pdfs = pdf.((cat_preference,), black_rat_center_mass)
black_rat_trap_rate
black_rat_catch_rates = black_rat_pdfs .* black_rat_trap_rate
black_rat_normalised_rates = black_rat_catch_rates ./ sum(black_rat_catch_rates)
black_rat_normalised_mass_yield = sum(black_rat_normalised_rates .* black_rat_mass_distribution .* mean_black_rat_mass)
black_rat_mass_yield = sum(black_rat_catch_rates .* black_rat_mass_distribution)
black_rat_predation_rate = sum(black_rat_catch_rates)

mouse_pdfs = pdf.((cat_preference,), mouse_center_mass)
mouse_trap_rate
mouse_catch_rates = mouse_pdfs .* mouse_trap_rate
mouse_normalised_rates = mouse_catch_rates ./ sum(mouse_catch_rates)
mouse_normalised_mass_yield = sum(mouse_normalised_rates .* mouse_mass_distribution .* mean_mouse_mass)
mouse_mass_yield = sum(mouse_catch_rates .* mouse_mass_distribution)
mouse_predation_rate = sum(mouse_catch_rates)

norway_rat_predation_rate / black_rat_predation_rate
norway_rat_predation_rate / mouse_predation_rate
black_rat_predation_rate / mouse_predation_rate
means = (mean_norway_rat_mass, mean_black_rat_mass, mean_mouse_mass)

predation_rates = norway_rat_predation_rate, black_rat_predation_rate, mouse_predation_rate
Ds = predation_rates ./ predation_rates[1]


# norway_rat : black_rat : mouse
# Roughly ~ 1:3:10
#
# rodents are ~80% of cat diet and they do not switch to natives much (Harper thesis etc)

fig = Figure(; size=(800, 600));
ax1 = Axis(fig[1, 1]; ylabel="Fraction trapped")
ax2 = Axis(fig[2, 1]; ylabel="Fraction killed")
ax3 = Axis(fig[3, 1]; xlabel="Prey size", ylabel="Reported means")
ax4 = Axis(fig[4, 1]; ylabel="Probability")
axs = ax1, ax2, ax3, ax4
hidexdecorations!.(axs[1:3]; grid=false)
hidespines!.(axs)
linkxaxes!(axs...)
alpha = 0.7
b1 = Makie.lines!(ax1, norway_rat_center_mass, norway_rat_trap_rate;
    color=(colors[1], alpha), label=rodent_labels[1] * " trapped sizes",
)
b2 = Makie.lines!(ax1, black_rat_center_mass, black_rat_trap_rate;
    color=(colors[2], alpha), label=rodent_labels[2] * " trapped sizes",
)
b3 = Makie.lines!(ax1, mouse_center_mass, mouse_trap_rate;
    color=(colors[3], alpha), label=rodent_labels[3] * " trapped sizes",
)
b4 = Makie.lines!(ax2, center_rat_size, killed_rats ./ sum(killed_rats);
    color=(colors[1], alpha), label="Norway rat killed sizes",
)
d = Makie.density!(ax3, cat_mean_prey_sizes; color=:grey, label="Literature mean prey sizes")
l = Makie.plot!(ax4, cat_preference; color=:black, label="Preference model")
axislegend.(axs; position=:rc)
Makie.xlims!.(axs, ((0, 600),))
display(fig)
save("images/cat_rodent_predation.png", fig)

# Makie.lines!(ax3, center_rat_size, killed_rats_glass ./ sum(killed_rats_glass))
# Makie.barplot!(ax3, center_rat_size, trapped_rats_childs ./ sum(trapped_rats_childs))
# Makie.lines!(ax3, center_rat_25_size, killed_rats_25_childs ./ sum(killed_rats_25_childs) .* 4)
# display(fig)


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

predator_fr_half_saturation = NamedVector(norway_rat=3, black_rat=5, mouse=8)
ks = Tuple(rodent_carrycap .* u"ha")
Ns = Tuple(rodent_carrycap .* u"ha") ./ 5
Ds = Tuple(predator_fr_half_saturation)
βs = D1 ./ Ds

c1 = cs[1]
c2 = cs[2]
c3 = cs[3]
k1 = ks[1]
k2 = ks[2]
k3 = ks[3]
r1 = rt[1]
r2 = rt[2]
r3 = rt[3]
N1 = Ns[1]
N2 = Ns[2]
N3 = Ns[3]
D1 = Ds[1]
D2 = Ds[2]
D3 = Ds[3]
β1 = βs[1]
β2 = βs[2]
β3 = βs[3]

RNV = NamedVector{(:norway_rat, :black_rat, :mouse),3}
t = 1u"yr"/12

N1, N2 = Ns
k1, k2 = ks
N1_ref = max(0.0, N1 + (r1 * N1 * (1 - (N1 + α12 * N2) / k1)) * 1u"yr" / 12  - ((c1 * P * β1 * N1 / (D1 + β1 * N1 + β2 * N2))) * 1u"yr"/12)
N2_ref = max(0.0, N2 + (r2 * N2 * (1 - (N2 + α21 * N1) / k2)) * 1u"yr" / 12  - ((c2 * P * β2 * N2 / (D1 + β1 * N1 + β2 * N2))) * 1u"yr"/12)

cs = Tuple(max_yield_fraction) ./ t
cat_pr0ference = LogNormal(log(50), 0.8)
ys = max_yield_fraction ./ t
ks = (25, 40, 100)
α12 = 0.5
α13 = 0.1
α21 = 2.0
α23 = 0.2
α31 = 4.0
α32 = 3.0
# αs = (x -> 1.0).(αs)
αs = @SMatrix [0.0 α12 α13; α21 0.0 α13; α31 α32 0.0]

# cs = individuals_per_cat
# Use relative predation rates for D ratios
# Ds = predation_rates ./ predation_rates[1]

hunted_rodent_mass = NamedVector(norway_rat=110, black_rat=90, mouse=15) .* u"g"
assimilated_energy_per_individual = hunted_rodent_mass .* rodent_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual
mean_prey_size = 60u"g" 
mean_prey_n = cat_energy_intake / (mean_prey_size * rodent_energy_content * assimilation_efficiency)

P = 0.00001
cs = individuals_per_cat
ys = max_yield_fraction ./ t
v = eachrow(pred_df)[1].rmax * u"yr^-1"
Ds = (4, 4, 4)
q = eachrow(pred_df)[1].carrycap
e = cat_energy_intake
Es = assimilated_energy_per_individual
Ncrit = 4 / u"d"
d_high = 0.2 / t
P = 0.01
Ns1 = 2 .* Ns
Ds = (1.0, 1.0, 1.0)
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

