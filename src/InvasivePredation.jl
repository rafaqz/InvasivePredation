module InvasivePredation

using Distributions
using NLopt
using Optimization
using OptimizationOptimJL
using OptimizationNLopt
using Unitful

export interaction_matrix, plot_densities, plot_interactions

export hanski_pred, hanski_predation, hanski_growth, hanski_multi

export find_predation_preference, optimize_predation_preference, 
    calculate_predation_rate, optimize_predation_rates_from_literature

export optimise_hunting

const basepath = realpath(joinpath(@__DIR__, ".."))

include("interactions.jl")
include("hanski.jl")
include("preference.jl")
include("yield.jl")

end
