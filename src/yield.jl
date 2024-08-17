

function rodent_func(x, p)
    (; k, rt, nsteps, years, seasonality, fixed_take, replicates) = p
    replicates = map(1:replicates) do i
        N = k
        # In a non-spatial model have to leave some population or it cant recover.
        minleft = oneunit(N) * 0.01
        total_taken = zero(N)

        # Burn in
        for _ in 1:5
            for s in 1:nsteps
                # @show N
                k_season = _seasonal_k(k, seasonality, nsteps, s)
                yield = rand(LogNormal(log(x[1]), p.std))
                N = (N .* k_season) ./ (N .+ (k_season .- N) .* exp.(.-(rt / (nsteps * u"yr^-1"))))
                # println()
                # @show k_season N
                # Cant take the whole population
                taken = if fixed_take
                    max(-N + minleft, k * -yield)
                else
                    max(-N + minleft, N * -yield)
                end
                N += taken
            end
        end
        total_taken = zero(N)
        for _ in 1:years
            for s in 1:nsteps
                # @show N
                k_season = _seasonal_k(k, seasonality, nsteps, s)
                yield = rand(LogNormal(log(x[1]), p.std))
                N = (N .* k_season) ./ (N .+ (k_season .- N) .* exp.(.-(rt / (nsteps * u"yr^-1"))))
                # println()
                # @show k_season N
                # Cant take the whole population
                taken = if fixed_take
                    max(-N + minleft, k * -yield)
                else
                    max(-N + minleft, N * -yield)
                end
                @assert taken <= zero(taken) "taken: $taken N: $N yield: $yield"
                @assert taken >= -N "taken: $taken N: $N yield: $yield"
                N += taken
                total_taken += taken
            end
        end
        return (taken=total_taken / years, final_pop=N)
    end
    return (taken=mean(first, replicates), final_pop=mean(last, replicates))
end
rodent_take_unitless(x, p) = ustrip(rodent_take(x, p))
rodent_take(x, p) = rodent_func(x, p)[1]
rodent_pop(x, p) = rodent_func(x, p)[2]

_seasonal_k(k, seasonality, nsteps, s) = k + k * seasonality * sin(2π * s / nsteps)

function optimise_hunting(rodent_params)
    p = first(rodent_params)
    optimal_takes = map(rodent_params) do p
        x0 = [0.1]
        prob = OptimizationProblem(rodent_take_unitless, x0, p; lb=[0.0], ub=[1.0])
        sol = solve(prob, SAMIN(); maxiters=10000)
        sol.u
    end
    optimal_pops = map(rodent_pop, optimal_takes, rodent_params)
    optimal_caught = map(rodent_take, optimal_takes, rodent_params) ./ u"yr" .* -1
    max_supported_cats = uconvert.(u"km^-2", optimal_caught ./ p.individuals_per_cat .* p.fraction_eaten)

    takes = NamedTuple{map(k -> Symbol(:taken_, k), propertynames(optimal_takes))}(values(optimal_takes))
    cats = NamedTuple{map(k -> Symbol(:cat_, k), propertynames(max_supported_cats))}(values(max_supported_cats))

    return (; std=p.std, seasonality=p.seasonality, map(first, takes)..., cats..., NamedTuple(optimal_pops)...)
end
