

function rodent_func(x, p)
    (; k, rt, nsteps, years, seasonality, fixed_take, replicates) = p
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
    (; k, rt, nsteps, fixed_take) = p
    k_season = _seasonal_k(k, seasonality, nsteps, s)
    N = (N .* k_season) ./ (N .+ (k_season .- N) .* exp.(.-(rt / (nsteps * u"yr^-1"))))
    yield = max(zero(x), rand(Normal(x, p.std * x)))
    # Cant take the whole population
    taken = min(N * yield, N - minleft)
    N -= taken
    # @assert taken <= zero(taken) "taken less than zero: $taken N: $N yield: $yield"
    # @assert taken >= -N "taken more than N: $taken N: $N yield: $yield"
    return N, taken
end

_seasonal_k(k, seasonality, nsteps, s) = k + k * seasonality * sin(2π * s / nsteps)

function optimise_hunting(rodent_params)
    optimal_fraction = map(rodent_params) do p
        x0 = [0.01]
        prob = OptimizationProblem(rodent_loss, x0, p; lb=[0.0], ub=[1.0])
        sol = solve(prob, SAMIN(); maxiters=10000)
        sol.u[1]
    end

    stats = map(optimal_fraction, rodent_params) do of, p
        mean(p.replicates) do _
            x = rodent_func(of, p)
            frac = of
            take = x.taken
            pop = x.mean_pop
            cats = uconvert(u"km^-2", pop * frac / p.individuals_per_cat / (1/12)u"yr")
            [frac, take, pop, cats]
        end |> NamedTuple{(:frac, :take, :pop, :cats)}
    end

    fracs = NamedTuple{map(k -> Symbol(:frac_, k), propertynames(rodent_params))}(map(x -> x.frac, stats))
    takes = NamedTuple{map(k -> Symbol(:taken_, k), propertynames(rodent_params))}(map(x -> x.take, stats))
    cats = NamedTuple{map(k -> Symbol(:cat_, k), propertynames(rodent_params))}(map(x -> x.cats, stats))
    pops = NamedTuple(map(x -> x.pop, stats))

    return (; std=rodent_params.mouse.std, seasonality=rodent_params.mouse.seasonality, fracs..., cats..., pops..., takes...)
end
