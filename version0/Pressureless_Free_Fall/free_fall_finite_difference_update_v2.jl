




function update_free_fall!(ρ, v, r, g, Δr, Δt)
    N = length(ρ)
    ρ_new = similar(ρ)
    v_new = similar(v)

    # Precompute conserved variables
    ρv   = ρ .* v
    ρvr2 = ρ .* v .* r.^2
    ρv2r2 = ρ .* v.^2 .* r.^2

    # Center (i = 1)
    i = 1
    ρ_new[i] = ρ[i] - (Δt / Δr) * (1 / r[i+1]^2) * (ρvr2[i+1] - ρvr2[i])
    v_num = ρv[i] - (Δt / Δr) * (1 / r[i+1]^2) * (ρv2r2[i+1] - ρv2r2[i]) - ρ[i] * g[i]
    v_new[i] = v_num / ρ_new[i]

    # Interior (2 ≤ i ≤ N-1)
    for i in 2:N-1
        ρ_new[i] = ρ[i] - (Δt / (2Δr)) * (1 / r[i]^2) * (ρvr2[i+1] - ρvr2[i-1])
        if ρ_new[i] < 0.0
            ρ_new[i] = 0.0
        end
        v_num = ρv[i] - (Δt / (2Δr)) * (1 / r[i]^2) * (ρv2r2[i+1] - ρv2r2[i-1]) - ρ[i] * g[i]
        v_new[i] = v_num / ρ_new[i]
    end

    # Surface (i = N)
    i = N
    ρ_new[i] = ρ[i] - (Δt / Δr) * (1 / r[i]^2) * (ρvr2[i] - ρvr2[i-1])
    v_num = ρv[i] - (Δt / Δr) * (1 / r[i]^2) * (ρv2r2[i] - ρv2r2[i-1]) - ρ[i] * g[i]
    v_new[i] = v_num / ρ_new[i]

    # Update arrays in-place
    ρ .= ρ_new
    v .= v_new
    return ρ, v
end
