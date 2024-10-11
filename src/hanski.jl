function hanski_sim(model)
    (; P, Ns, P_timeline, Ns_timeline, ks, α, tspan) = model
    P1 = P
    Ns1 = Ns
    @inbounds for i in tspan
        Ns2 = hanski_prey_timestep(P1, Ns1, ks, α, model)::typeof(Ns)
        P1 = hanski_predator_timestep(P1, Ns1, model)::typeof(P)
        Ns1 = Ns2
        P_timeline[i] = P1
        Ns_timeline[i] = Ns1
    end
    return P_timeline, Ns_timeline
end

# Multi-prey weighted predation
@inline predation(N::Number, β::Number, c::Number, P::Number, D_Nβs::Number) =
    c * P * β * N / D_Nβs

# Simple logistic population growth model
@inline growth(N::Number, k::Number, r::Number, Ns_x, α_x) =
    r * N * (1 - (N + sum(α_x .* Ns_x)) / k)

# TODO generate this
const OTHER_INDS = ((2, 3), (1, 3), (1, 2))

# N-prey predation model
@inline function hanski_prey_timestep(P, Ns::NamedVector{K}, ks, α, model) where K
    (; Ds, cs, rs, t) = model
    @inbounds βs = map(d -> Ds[1] / d, Ds)
    Nβs = sum(βs .* Ns)
    @inbounds D_Nβs = Ds[1] + Nβs
    is = NamedVector{K}(ntuple(identity, Val{length(K)}()))
    @inbounds map(is) do i
        other_inds = OTHER_INDS[i]
        N = Ns[i]
        α_x = α[i]
        Ns_x = map(j -> Ns[j], other_inds)
        g = growth(Ns[i], ks[i], rs[i], Ns_x, α_x)
        p = predation(Ns[i], βs[i], cs[i], P, D_Nβs)
        max(zero(N), N + (g - p) * t)
    end
end

# Predator growth rate is independent from hunting
# the breeding threshold is like perception of
# excess rather than current intake
@inline function hanski_predator_timestep(P, Ns, model)
    (; v, e, d_high, ys, Es, Ds, t) = model
    βs = Ds[1] ./ Ds
    Nβs = sum(Ns .* βs)
    # Calculate carrycap `q` from yield and energy requirement
    q = convert(typeof(P), sum(Ns .* ys .* Es) / e)
    Preproduction = 1.5P # ??
    # If prey are above the breeding threshold
    P1 = if q > Preproduction # Use supportable population as the threshold rather than Ncrit
        max(zero(P), P + v * P * (1 - q * P / Nβs) * t)
    else
        max(zero(P), P + -d_high * P * t)
    end
    P1 + model.stochasticity * (rand() - 0.5) * P1
end
