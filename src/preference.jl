
function find_mass_distribution(x, p)
    μ, σ = x
    mass_distibution = LogNormal(log(μ), σ)
    halfstep = step(p.bin_center_mass) / 2
    predicted_rate = cdf.(mass_distibution, p.bin_center_mass .+ halfstep) .-
        cdf.(mass_distibution, p.bin_center_mass .- halfstep)
    # Get the rest of the distribution
    empty_upper = cdf(mass_distibution, Inf) - cdf(mass_distibution, last(p.bin_center_mass) + 1halfstep)
    residuals = p.trap_rate .- predicted_rate
    # Return the sum of squared residuals including the empty top end
    return sum(residuals .^ 2) + empty_upper .^ 2
end

function find_predation_preference(x, p)
    μ, σ = x
    cat_preference = LogNormal(log(μ), σ);
    predation_pds = pdf.((cat_preference,), p.bin_center_mass) .* p.trapped_rats
    killed_means = p.killed_rats ./ sum(p.killed_rats)
    predicted_means = predation_pds ./ sum(predation_pds)

    return sum((killed_means .- predicted_means) .^ 2)
end

function calculate_predation_rate(predator_preference, bin_center_mass, trap_rate)
    # Calculate bin mass statistics
    bin_mean_mass = trap_rate .* bin_center_mass
    mean_mass = sum(bin_mean_mass)
    mass_distribution = bin_mean_mass ./ mean_mass

    # Calculate predation probability densities for mass bins
    # based on the mass preference distribution of the predator
    predation_pds = pdf.((predator_preference,), bin_center_mass)

    # We assume the trap rate is similar to 
    # the encounter rate experience by predators
    predation_rates = predation_pds .* trap_rate
    normalised_rates = predation_rates ./ sum(trap_rate)
    normalised_mass_yield = sum(normalised_rates .* mass_distribution .* mean_mass)
    mass_yields = sum(predation_rates .* mass_distribution)
    predation_rate = sum(predation_rates)
    return (; bin_center_mass, mean_mass, mass_distribution, mass_yields, predation_pds, predation_rates, predation_rate, trap_rate)
end

function optimize_mass_distribution(p)
    x0 = [1.0, 1.0]
    find_mass_distribution(x0, p)
    prob = OptimizationProblem(find_mass_distribution, x0, p;
        lb=[0.0, 0.0],
        ub=[1000.0, 1000.0]
    )
    sol = solve(prob, NLopt.LN_NELDERMEAD())
    μ, σ = sol.u
    return LogNormal(log(μ), σ)
end

function optimize_predation_preference(p)
    x0 = [1.0, 1.0]
    find_predation_preference(x0, p)
    prob = OptimizationProblem(find_predation_preference, x0, p;
        lb=[0.0, 0.0],
        ub=[500.0, 10.0]
    )
    sol = solve(prob, NLopt.LN_NELDERMEAD())
    μ, σ = sol.u
    return LogNormal(log(μ), σ)
end

function fit_distributions_to_literature()
    # Bins
    bin_max_100 = 100:100:700 # bin max - Trapped rats binned to 100g
    bin_center_100 = 50:100:650 # bin centers
    bin_center_25 = 12.5:25:187.5 # Killed rats binned to 25g

    # Counts
    killed_rats_25_childs = [0, 8, 8, 3, 1, 1, 0, 1]
    killed_rats_childs = [19, 4, 0, 0, 0, 0, 0]
    trapped_rats_childs = [5, 5, 9, 29, 44, 14, 3]

    # sum(center_rat_25_size .* killed_rats_25_childs) ./ sum(killed_rats_25_childs)

    # Data form Glass 2009

    # Counts
    trapped_rats_glass = [18, 38, 49, 167, 175, 98, 0]
    killed_rats_glass = [15, 13, 4, 2, 1, 0, 0]

    mean_glass = sum(bin_center_100 .* killed_rats_glass) ./ sum(killed_rats_glass)
    n_glass = sum(killed_rats_glass)

    mean_childs = sum(bin_center_25 .* killed_rats_25_childs) ./ sum(killed_rats_25_childs)
    n_childs = sum(killed_rats_25_childs)

    n_total = n_glass + n_childs
    frac_glass = n_glass / n_total
    frac_childs = n_childs / n_total
    mean_caught_norway_rat_mass = mean_glass .* frac_glass + mean_childs .* frac_childs

    norway_rat_studies = (
        bin_center_25,
        bin_center_100,
        glass = (;
            trapped=trapped_rats_glass,
            killed=killed_rats_glass,
            mean=mean_glass,
            n=n_glass,
        ),
        childs = (;
            trapped=trapped_rats_childs,
            killed=killed_rats_childs,
            killed_25=killed_rats_25_childs,
            mean=mean_childs,
            n=n_childs,
        ),
    )

    trapped_rats = trapped_rats_glass .+ trapped_rats_childs
    killed_rats = killed_rats_glass .+ killed_rats_childs

    # Model parameters
    norway_rat_params = (; trapped_rats, killed_rats, bin_center_mass=bin_center_100)

    # Optimise a LogNormal distribution
    cat_mass_preference = optimize_predation_preference(norway_rat_params)

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

    # Wilson et al 2007
    black_rat_bin_center_mass = (40:40:240) .- 20
    black_rat_trap_rate = [0.08, 0.20, 0.26, 0.31, 0.14, 0.01]

    # Glass 1988
    norway_rat_bin_center_mass = (100:100:500) .- 50
    male_norway_rat_trap_rate = [0.07, 0.3, 0.43, 0.16, 0.04]
    female_norway_rat_trap_rate = [0.14, 0.33, 0.34, 0.19, 0.0]
    norway_rat_trap_rate = (male_norway_rat_trap_rate .+ female_norway_rat_trap_rate) ./ 2

    # Estimated from shifted rat samples
    mouse_bin_center_mass = 5:10:55
    mouse_trap_rate = [0.08, 0.20, 0.26, 0.31, 0.14, 0.01]

    trap_rates = (;
        norway_rat=norway_rat_trap_rate,
        black_rat=black_rat_trap_rate,
        mouse=mouse_trap_rate,
    )
    bin_center_masses = (;
        norway_rat=norway_rat_bin_center_mass,
        black_rat=black_rat_bin_center_mass,
        mouse=mouse_bin_center_mass,
    )

    rodent_mass_distributions = map(trap_rates, bin_center_masses) do trap_rate, bin_center_mass
        optimize_mass_distribution((; trap_rate, bin_center_mass))
    end

    rodent_stats =  map(bin_center_masses, trap_rates) do bcm, tr
        calculate_predation_rate(cat_mass_preference, bcm, tr)
    end

    return (; cat_mass_preference, rodent_mass_distributions, rodent_stats, norway_rat_params, norway_rat_studies)
end
