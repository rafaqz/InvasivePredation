function hanski_predation(N::Number, β::Number, c::Number, P::Number, D_Nβs::Number)
    # Multi-prey weighted predation
    c * P * β * N / D_Nβs
end

function hanski_growth(N::Number, k::Number, r::Number, Ns_x, αs_x)
    r * N * (1 - (N + sum(αs_x .* Ns_x)) / k)
end

function hanski_multi(P::Number, Ns::NTuple{I}, Ds, Es, ys, αs, ks, cs, rs, d_high, t) where I
    βs = Ds[1] ./ Ds
    Nβs = sum(βs .* Ns)
    D_Nβs = Ds[1] + Nβs
    is = ntuple(identity, Val{I}())
    return map(is) do i
        N = Ns[i]
        others = Tuple([j for j in is if j != i])
        αs_x = map(j -> αs[i, j], others)
        Ns_x = map(j -> Ns[j], others)
        growth = hanski_growth(Ns[i], ks[i], rs[i], Ns_x, αs_x)
        predated = hanski_predation(Ns[i], βs[i], cs[i], P, D_Nβs)
        N1 = max(zero(N), N + (growth - predated) * t)
        return N1
    end
end

function hanski_pred(P::Number, v::Number, e::Number, d_high::Number, Ns, ys, Es, Ds, cs, t)
    βs = Ds[1] ./ Ds
    Nβs = sum(Ns .* βs)
    q = convert(typeof(P), sum(Ns .* ys .* Es) / e)
    Preproduction = 1.5P # ??
    # If prey are above the breeding threshold
    P1 = if q > Preproduction # Use supportable population as the threshold rather than Ncrit
        max(zero(P), P + v * P * (1 - q * P / Nβs) * t)
    else
        max(zero(P), P + -d_high * P * t)
    end
    P1
end
