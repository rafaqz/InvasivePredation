using Revise
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
using CairoMakie
CairoMakie.activate!()

using InvasivePredation

basepath = InvasivePredation.basepath

s = InvasivePredation.load_settings()
(; stochastic_rates, max_yield_fraction) = InvasivePredation.get_max_yield_fraction()

# Stochasticity reduces yield
fig = let
    fig = Figure(; size=(500, 500));
    labels = ["Optimal fraction\nhunted per month", "Supported cats\n", "Mean rodents\n"]
    yl = ((0, 0.5), (0, 2.0), (0, 40))
    x = stochastic_rates.std
    xticks = 0:0.4:4.0
    linewidth = 2
    ax_kw = (; xticks, spinewidth=2, xgridwidth=2, ygridwidth=2)
    ax_yticks = (0.1:0.1:0.4, 0.5:0.5:1.5, 10:10:30)
    axs = map(enumerate(labels)) do (i, ylabel)
        yticks = ax_yticks[i]
        ax = if i == 3
            Axis(fig[i, 1];
                xlabel="Scalar of standard deviation\nof hunting stochasticity", 
                yticks,
                ylabel, ax_kw...
            )
        else
            Axis(fig[i, 1]; 
                yticks,
                xticksvisible=false, 
                xticklabelsvisible=false,
                ylabel, ax_kw...
            )
        end
        label!(fig, ('A':'E')[i], i, 0)
        hidexdecorations!(ax; label=false, ticks=false, ticklabels=false)
        hideydecorations!(ax; label=false, ticks=false, ticklabels=false)
        ylims!(ax, yl[i])
        ax
    end
    colsize!(fig.layout, 0,  Relative(0.05))
    linkxaxes!(axs...)
    rodents = stochastic_rates[:, rodent.names]
    fracs = stochastic_rates[:, map(x -> Symbol(:frac_, x), rodent.names)]
    cats = stochastic_rates[:, map(x -> Symbol(:cat_, x), rodent.names)]
    rodent_plots = map(1:3) do i
        color = InvasivePredation.rodent_colors[i]
        label = s.rodent.labels[i]
        p1 = Makie.lines!(axs[1], x, fracs[:, i];
            label, color, linewidth
        )
        p2 = Makie.lines!(axs[2], x, cats[:, i];
            label, color, linewidth
        )
        p3 = Makie.lines!(axs[3], x, rodents[:, i];
            label, color, linewidth
        )
        (p1, p2, p3)
    end
    axislegend(axs[1])
    fig
end

save(joinpath(basepath, "images/cat_hunting_stochasticty.png"), fig)
