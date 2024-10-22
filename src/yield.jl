const MONTH = (1/12)u"yr"

function get_max_yield_fraction()
    (; cat, rodent, yield_metaparams) = InvasivePredation.load_settings()
    y = yield_metaparams

    (; cat_mass_preference, rodent_stats, rodent_mass_distributions, norway_rat_params, norway_rat_studies) = 
        fit_distributions_to_literature()

    (; individuals_per_cat) = InvasivePredation.get_cat_energetics(cat, rodent, rodent_stats)
    # Note: we get the same yield regardless of k, so this is general to any k
    st = 0.0
    stochastic_rates = map((0.0:0.4:4.0)) do st
        rodent_params = map(rodent.rmax, rodent.carrycap, individuals_per_cat) do rt, k, ipc
            (; k, rt, std=st, metaparams=yield_metaparams, individuals_per_cat=ipc)
        end
        stats = map(optimise_hunting, rodent_params)
        fracs = NamedTuple{map(k -> Symbol(:frac_, k), propertynames(rodent_params))}(map(x -> x.frac, stats))
        takes = NamedTuple{map(k -> Symbol(:taken_, k), propertynames(rodent_params))}(map(x -> x.take, stats))
        cats = NamedTuple{map(k -> Symbol(:cat_, k), propertynames(rodent_params))}(map(x -> x.cats, stats))
        pops = NamedTuple(map(x -> x.pop, stats))

        return (; std=rodent_params.mouse.std, fracs..., cats..., pops..., takes...)
    end |> DataFrame
    max_yield_fraction = get_yield_fraction(stochastic_rates, rodent, 1)
    return (; stochastic_rates, max_yield_fraction)
end

function get_yield_fraction(df, rodent, i)
    NamedVector{Tuple(rodent.names)}(
        collect(df[1, map(x -> Symbol(:frac_, x), rodent.names)])
    ) * 12 * u"yr^-1"
end

function get_cat_energetics(cat, rodent, rodent_stats)
    mean_preferred_mass = NamedVector(map(r -> r.mean_predation_mass, rodent_stats))
    assimilated_energy = mean_preferred_mass .* rodent.energy_content .* cat.assimilation_efficiency .* cat.fraction_eaten
    individuals_per_cat = cat.energy_intake ./ assimilated_energy
    return (; assimilated_energy, individuals_per_cat, mean_preferred_mass)
end

function optimise_hunting(p)
    x0 = [0.01]
    prob = OptimizationProblem(rodent_loss, x0, p; lb=[0.0], ub=[1.0])
    sol = solve(prob, SAMIN(); maxiters=10000)
    frac = sol.u[1]
    x = rodent_func(frac, p)
    take = x.taken
    pop = x.mean_pop
    cats = uconvert(u"km^-2", pop * frac / p.individuals_per_cat / MONTH)
    return (; frac, take, pop, cats)
end

function rodent_func(x, p)
    (; k, rt, metaparams) = p
    (; nsteps, years, seasonality, replicates) = metaparams
    reps = map(1:replicates) do i
        N = 2k
        # In a non-spatial model have to leave some population or it cant recover.
        minleft = oneunit(N) * 0.01
        nyears = length(oneunit(years):oneunit(years):years)

        # Burn in
        for _ in 1:nyears
            for s in 1:nsteps
                N, taken = _update_pop(seasonality, s, x[1], minleft, p, N)
            end
        end

        sumN = zero(N)
        total_taken = zero(N)
        for _ in 1:nyears
            for s in 1:nsteps
                sumN += N
                N, taken = _update_pop(seasonality, s, x[1], minleft, p, N)
                total_taken += taken
            end
        end
        (taken=total_taken / years, final_pop=N, mean_pop=sumN / (nsteps * nyears))
    end
    return (taken=mean(r -> r.taken, reps), final_pop=mean(r -> r.final_pop, reps), mean_pop=mean(r -> r.mean_pop, reps))
end
rodent_loss(x, p) = -ustrip(rodent_func(x, p)[:taken])

function _update_pop(seasonality, s, x, minleft, p, N)
    (; k, rt, metaparams) = p
    (; nsteps) = metaparams
    k_season = _seasonal_k(k, seasonality, nsteps, s)
    N = (N .* k_season) ./ (N .+ (k_season .- N) .* exp.(.-rt))
    yield = max(zero(x), rand(Normal(x, p.std * x)))
    # Cant take the whole population
    taken = min(N * yield, N - minleft)
    @assert taken >= zero(taken) "taken less than zero: $taken N: $N yield: $yield"
    @assert taken <= N "taken more than N: $taken N: $N yield: $yield"
    N -= taken
    return N, taken
end

_seasonal_k(k, seasonality, nsteps, s) = k + k * seasonality * sin(2π * s / nsteps)
