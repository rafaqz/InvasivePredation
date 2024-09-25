module InvasivePredation

using Distributions
using NLopt
using Optimization
using OptimizationOptimJL
using OptimizationNLopt
using Unitful
using LandscapeChange

export interaction_matrix, plot_densities, plot_interactions

export hanski_predator_timestep, hanski_prey_timestep

export find_predation_preference, optimize_predation_preference, 
    calculate_predation_rate, fit_distributions_to_literature

export optimise_hunting

const basepath = realpath(joinpath(@__DIR__, ".."))

include("interactions.jl")
include("hanski.jl")
include("preference.jl")
include("yield.jl")

end
