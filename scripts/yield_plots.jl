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

# Use mean from actual sizes taken for Norway rats, and guess the others
# mean_rodent_mass = R(rodent_df.mass) .* u"g"
#
#aTODO calculate these means using the preference model
#
(; cat_mass_preference, rodent_stats, rodent_mass_distributions, norway_rat_params, norway_rat_studies) = 
    fit_distributions_to_literature()

# Use the mean selected mass calculated earlier
mean_preferred_rodent_mass = NamedVector(map(r -> r.mean_predation_mass, rodent_stats))
assimilated_energy_per_individual = mean_preferred_rodent_mass .* rodent_energy_content .* assimilation_efficiency .* fraction_eaten
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual


function get_max_yield_fraction(p)
    stochastic_rates = map((0.0:0.4:4.0)) do st
        rodent_params = map(rodent_rmax, rodent_carrycap, individuals_per_cat) do rt, k, ipc
            (; k, rt, replicates, fixed_take=false, nsteps, years, std=st, seasonality, fraction_eaten, individuals_per_cat=ipc)
        end
        optimise_hunting(rodent_params)
    end |> DataFrame
    max_yield_fraction = NamedVector{Tuple(rodent_names)}(
        collect(stochastic_rates[1, map(x -> Symbol(:frac_, x), rodent_names)])
    )
end

# k does not matter - we get the same yield rate regardless of k
# so the model is general to any k

# Stochasticity reduces yield
fig = let
    fig = Figure(; title="Effects of monthly stochasticity on optimal hunting yields", size=(600, 800));
    labels = ["Optimal fraction hunted per month", "Supported cats per km^2", "Mean rodents per hectare"]
    yl = ((0, 0.2), (0, 1.5), (0, 40))
    x = stochastic_rates.std
    xticks = 0:0.4:4.0
    axs = map(enumerate(labels)) do (i, ylabel)
        ax = if i == 3
            Axis(fig[i, 1];
                xlabel="Scalar of standard deviation of hunting stochasticity", 
                # xlabel="seasonality", 
                ylabel, 
                xticks,
            )
        else
            Axis(fig[i, 1]; 
                ylabel, 
                xticks,
                xticksvisible=false, 
                xticklabelsvisible=false
            )
        end
        ylims!(ax, yl[i])
        # hidespines!(ax)
        ax
    end
    linkxaxes!(axs...)
    rodents = stochastic_rates[:, rodent_names]
    fracs = stochastic_rates[:, map(x -> Symbol(:frac_, x), rodent_names)]
    cats = stochastic_rates[:, map(x -> Symbol(:cat_, x), rodent_names)]
    rodent_plots = map(1:3) do i
        color = colors[i]
        label = rodent_labels[i]
        p1 = Makie.lines!(axs[1], x, fracs[:, i];
            label, color
        )
        p2 = Makie.lines!(axs[2], x, cats[:, i];
            label, color
        )
        p3 = Makie.lines!(axs[3], x, rodents[:, i];
            label, color
        )
        (p1, p2, p3)
    end
    axislegend(axs[1])
    fig
end

save(joinpath(basepath, "images/cat_hunting_stochasticty.png"), fig)
