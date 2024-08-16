# convenience functions
# picks n random numbers between low and high, that are no more than a factor 2 apart
function _rejection_range(n, low=0, high=1, max_multiplier=2)
    while true
        cont = false
        r = rand(n) .* (high - low) .+ low
        for i in 1:n, j in 1:n
            if r[i] / r[j] > max_multiplier || r[j] / r[i] > max_multiplier
                cont = true
            end
        end
        cont || return r
    end
end

function interaction_matrix(bodymass;
    competition_to_predation_max = 0.1
)
    # Create interaction matrix
    cats = 1
    norway_rat = 2
    black_rat = 3
    mouse = 4
    pig = 5
    rodents = [black_rat, norway_rat, mouse]

    interactions = zeros(5, 5)

    # cat/rodent
    # relative rankings of the effect of cats on an animal, low to high
    cat_order = [norway_rat, black_rat, mouse]
    interactions[cats, rodents] = -sort(rand(3))[sortperm(cat_order)]

    interactions[rodents, cats] .= (-interactions[cats, i] * sqrt(bodymass[i]) / sqrt(bodymass[cats]) for i in rodents)

    # cat/pig
    # Small negative interactions
    smallest_interaction = minimum(abs.(interactions[cats, black_rat:mouse]))
    # Pigs and cats are randomly less than half of the
    # smallest interactions of the other species
    interactions[cats, pig] = -rand() * 0.5smallest_interaction
    interactions[pig, cats] = -rand() * 0.5smallest_interaction

    # rodents
    ninteractions = 6
    low = 0
    high = competition_to_predation_max * abs(interactions[cats, black_rat])
    max_multiplier = 2
    randvals = -_rejection_range(ninteractions, low, high, max_multiplier)
    for i in black_rat:mouse, j in black_rat:mouse
        if i != j
            interactions[i, j] = pop!(randvals) * sqrt(bodymass[i]) / sqrt(bodymass[j])
        end
    end

    # pig/rodent
    ninteractions = 6
    low = 0
    high = minimum(abs(x) for x in interactions if x!= 0)
    max_multiplier = 2
    randpigs = -_rejection_range(ninteractions, low, high, max_multiplier)
    for i in black_rat:mouse
        interactions[i, 5] = pop!(randpigs)
        interactions[5, i] = pop!(randpigs)
    end

    return interactions
end

# visualize a single interaction matrix
function plot_interactions(mat, nodes)
    vals = [mat[j, i] for i in axes(mat, 1), j in axes(mat, 2) if i != j]
    cols = ifelse.(vals .>= 0, :blue, :red)
    g = complete_digraph(5)
    f = Figure(size = (800, 800))
    graphplot(f[1,1], g, layout = Shell(), nlabels = nodes, edge_width = 3abs.(vals), edge_color = cols)

    return f
end

# visualize an aggregate from a vector of interaction matrices
function plot_densities(mat, nodes)
    mm = stack(mat)

    f = Figure(size = (900, 1200))
    for i in 1:5, j in 1:5
        if i != j
            col = sum(mm[i, j, :]) < 0 ? :red : :blue
            ax = Axis(f[i,j];
                title="$(nodes[i]) on $(nodes[j])",
                limits=(-1, 1, nothing, nothing)
            )
            GLMakie.hist!(ax, mm[i,j,:],
                bins=200, color=(col, 0.3), strokearound=true, strokewidth=1, strokecolor=col
            )
        end
    end
    f
end
