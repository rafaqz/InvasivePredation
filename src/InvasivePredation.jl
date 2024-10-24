module InvasivePredation

using Distributions
using NLopt
using Optimization
using OptimizationOptimJL
using OptimizationNLopt
using Unitful
using LandscapeChange
using ColorSchemes
using CSV
using DataFrames

export interaction_matrix, plot_densities, plot_interactions

export hanski_predator_timestep, hanski_prey_timestep, hanski_sim

export find_predation_preference, optimize_predation_preference, 
    calculate_predation_rate, fit_distributions_to_literature

export optimise_hunting

const basepath = realpath(joinpath(@__DIR__, ".."))

include("hanski.jl")
include("preference.jl")
include("yield.jl")
include("load_settings.jl")

cat_color = ColorSchemes.solar[0.7]
rodent_colors = ColorSchemes.tempo[[0.2, 0.55, 0.9]]


end
