




function update_density!(ρ::AbstractVector, v::AbstractVector, r::AbstractVector,
                         Δt::Float64, Δr::Float64)
    N = length(ρ)
    ρ_new = zeros(N)
    @assert length(v)==N && length(r)==N
    # Precompute mass‐flux term f[i] = ρ[i]*v[i]*r[i]^2
    f = similar(ρ)
    @inbounds for i in 1:N
        f[i] = ρ[i] * v[i] * r[i]^2
    end

    # Update interior points (i=2..N-1)
    @inbounds for i in 2:N-1
        ρ_new[i-1] = ρ[i-1] - (Δt/(2Δr)) * (1/r[i]^2) * (f[i+1] - f[i-1])
    end

    return ρ_new
end

function update_velocity!(
    v::AbstractVector,
    g::AbstractVector,
    rho::AbstractVector,
    P::AbstractVector,
    Δt::Float64,
    Δr::Float64
)
    N = length(v)
    v_new = zeros(N)

    @inbounds for i in 2:N-1
        if rho[i] <= 0.0
            v_new[i-1] = v[i-1] + g[i-1]*Δt - (v[i-1]*Δt)/2Δr * (v[i+1] - v[i-1])
        else
            v_new[i-1] = v[i-1] + g[i-1]*Δt - (v[i-1]*Δt)/2Δr * (v[i+1] - v[i-1]) + (Δt/(2 * rho[i-1] * Δr)) * (P[i+1] - P[i-1])
        end
    end

    return v_new
end

function update_energy!(
    E::AbstractVector,
    ρ::AbstractVector,
    v::AbstractVector,
    g::AbstractVector,
    P::AbstractVector,
    Δt::Float64,
    Δr::Float64
)
    N = length(E)
    @assert length(ρ)==N && length(v)==N && length(g)==N && length(P)==N

    S = similar(E)
    F = similar(E)

    @inbounds for i in 1:N
        S[i] = Δt * ρ[i] * v[i] * g[i]
        F[i] = (E[i] + P[i]) * v[i]
    end

    E_tmp = copy(E)  # snapshot of original E values

    @inbounds for i in 2:N-1
        E[i] = E_tmp[i] +
               S[i] -
               (Δt / (2 * Δr)) * (F[i+1] - F[i-1])
    end

    # Boundary update (simple forward difference)
    E[1] = E_tmp[1] + S[1] - (Δt / Δr) * (F[2] - F[1])
    E[end] = E_tmp[end] + S[end] - (Δt / Δr) * (F[end] - F[end-1])
    
    return E
end
