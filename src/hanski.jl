function hanski_sim(model;
    P_timeline=Vector{typeof(model.P0)}(undef, length(model.tspan)),
    Ns_timeline=Vector{typeof(model.Ns0)}(undef, length(model.tspan)),
)
    (; P0, Ns0, ks, α, tspan) = model
    P = P0
    Ns = Ns0
    @inbounds for (i, t) in enumerate(tspan)
        Ns1 = Ns_timeline[i] = hanski_prey_timestep(P, Ns, ks, α, model)::typeof(Ns0)
        P = P_timeline[i] = hanski_predator_timestep(P, Ns, model)::typeof(P0)
        Ns = Ns1
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
        convert(eltype(Ns), max(zero(N), N + (g - p * t)))
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
    q = convert(typeof(P), (sum(Ns .* ys .* Es)) / e)
    @show q Ns ys Es e
    println()
    Preproduction = 1.2P # ??
    # If prey are above the breeding threshold
    P1 = if q > Preproduction # Use supportable population as the threshold rather than Ncrit
        max(zero(P), P + v * P * (1 - q * P / oneunit(P) / Nβs))
    else
        max(zero(P), P + -d_high * P * t)
    end
    P1 + model.stochasticity * (rand() - 0.5) * P1
end

