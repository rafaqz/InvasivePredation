using CSV
using DataFrames
using Distributions
using LandscapeChange
using Optimization
using OptimizationOptimJL
using StaticArrays
using Unitful
using GLMakie
using Distributions
using NLopt
using OptimizationNLopt
# using CairoMakie
# CairoMakie.activate!()

using InvasivePredation

basepath = InvasivePredation.basepath

include(joinpath(basepath, "scripts/load_settings.jl"))

# Use mean from actual sizes taken for Norway rats, and guess the others
mean_rodent_mass = R(rodent_df.mass) .* u"g"
# TODO calculate these means using the preference model
mean_preferred_rodent_mass = NamedVector(norway_rat=110, black_rat=90, mouse=15) .* u"g"
assimilated_energy_per_individual = mean_rodent_mass .* rodent_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual

# nsteps = 365
# 12 steps (months) is a close enough approximation of continuous

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

nsteps = 12
years = 10
fraction_eaten = 0.72 # McGregor 2015
seasonality = 0.5
replicates = 1

stochastic_rates = map((0.0:0.2:1.5).^2) do std
    rodent_params = map(R(rodent_names), individuals_per_cat) do rodent, individuals_per_cat
        (; k=k[rodent], rt=rt[rodent], replicates=25, fixed_take=false, nsteps, years, std, seasonality, fraction_eaten, individuals_per_cat)
    end
    optimise_hunting(rodent_params)
end |> DataFrame
max_yield_fraction = Tuple(DataFrame(stochastic_rates)[1, 3:5])

# Why is this relationship slightly under 1
# max_yield_fraction ./ (takes[1] / ustrip(rt[1]))

# k does not matter - we get the same yield rate regardless of k
# so the model is general to any k

# Stochasticity reduces yield
fig = let
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
    x = stochastic_rates.std
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
    fig
end

save(joinpath(basepath, "images/cat_hunting_stochasticty.png"), fig)
