




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
    Δt::Float64,
    Δr::Float64
)
    N = length(v)
    v_new = zeros(N)

    @inbounds for i in 2:N-1
        v_new[i] = v[i] + g[i]*Δt - (v[i]*Δt/(2Δr)) * (v[i+1] - v[i-1])
    end
    v_new[1] = 0.0 # Boundary condition at the center
    v_new[N] = 0.0
    return v_new
end