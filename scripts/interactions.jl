
using Graphs, GLMakie, GraphMakie, NetworkLayout, InvertedIndices, LandscapeChange
using CSV, DataFrames, StaticArrays

# Metadata

pred_names = ["cat", "black_rat", "norway_rat", "mouse", "pig"]
# run(`libreoffice tables/animals.csv`)
pred_df = filter(CSV.read("tables/animals.csv", DataFrame)) do row
    row.name in pred_names
end

bodymass = map(identity, pred_df.mass)

# one matrix
ints = interaction_matrix(bodymass)
plot_interactions(ints, pred_names)

# pred_pops = NamedVector{Tuple(map(Symbol, pred_names)),5}(rand(5)) .* 10
# pred_interactions = SMatrix{5,5}(interaction_matrix(bodymass))
# pred_interactions * pred_pops

# multiple
mat = map(_-> interaction_matrix(bodymass), 1:1000)
plot_densities(mat, pred_names)
map(1:100) do _
    interaction_matrix(bodymass)
end
