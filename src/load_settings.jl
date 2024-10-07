function load_settings()
    animals_csv = joinpath(InvasivePredation.basepath, "tables/animals.csv")
    # run(`libreoffice $animals_csv`)
    pred_df = CSV.read(animals_csv, DataFrame)
    rodent_df = pred_df[2:4, :]
    rodent_names = map(Symbol, rodent_df.name)
    R = NamedVector{Tuple(rodent_names),3}
    cat_data = map(identity, NamedTuple(filter(r -> r.name == "cat", pred_df)[1, :]))

    return (;
        rodent=(
            energy_content=6.24u"kJ/g",
            carrycap=R(rodent_df.carrycap) .* u"ha^-1",
            rmax=R(rodent_df.rmax) .* u"yr^-1",
            colors=R([:red, :lightblue, :yellow]),
            names=rodent_names,
            labels=titlecase.(replace.(rodent_df.name, ("_" => " ",)))
        ),
        cat=(
            mean_prey_sizes=sort!(
                [11.0, 33.0, 17.9, 4.3, 15.5, 23.6, 24.0, 31.2, 26.5, 23.3, 9.3, 26.2, 29.0, 
                 21.2, 349.7, 3.64, 21.8, 35.9, 13.6, 8.2, 228.6, 249.0, 42.7, 16.5, 15.5, 
                 60.7, 32.8, 184.9, 34.1, 36.8, 7.4, 102.0, 241.2, 38.5, 13.3, 26.4, 7.6, 24.3, 
                 123.2, 39.2, 1.77, 76.7, 13.5, 18.3, 30.3, 34.4]
            ),
            rmax=cat_data.rmax,
            energy_intake=2131u"kJ/d",
            fraction_eaten=0.72, # McGregor 2015
            assimilation_efficiency=0.84,
        ),
        yield_metaparams=(;
            nsteps=12, # 12 steps (months) is a close enough approximation of continuous
            years=10u"yr",
            seasonality=0.0,
            replicates=10,
        )
    )
end
