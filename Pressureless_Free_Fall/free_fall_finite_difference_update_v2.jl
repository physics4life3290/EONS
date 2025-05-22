


# This method is not adequate #
# Needs more cowbell #
# This is a finite difference update for the free fall problem...its not good #

function update_freefall(
    ρ::AbstractVector{<:Real},
    v::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    g_acc::AbstractVector{<:Real},
    dt::Real,
    dr::Real
)
    N = length(ρ)
    @assert length(v)==N && length(r)==N && length(g_acc)==N "All input vectors must have the same length"

    # compute fluxes once
    ρ_flux = ρ .* v .* r .^ 2
    mom_flux = ρ_flux .* v

    # allocate outputs
    ρ_new = similar(ρ)
    v_new = similar(v)

    # interior + boundary
    for i in 1:N
        if i == 1
            # i=1: use one‐sided difference forward
            Δρ = (ρ_flux[2])
            ρ_new[1] = ρ[1] - (3 * dt / dr) * Δρ

            Δm = (mom_flux[2] - mom_flux[1]) / r[1]^2
            v_new[1] = (ρ[1]*v[1] - (dt/(2dr))*Δm - ρ[1]*g_acc[1]) / ρ_new[1]

        elseif i == N
            # outer boundary: zero out
            ρ_new[N] = 0.0
            v_new[N] = 0.0

        else
            # central difference for interior points
            Δρ = (ρ_flux[i+1] - ρ_flux[i-1]) / r[i]^2
            ρ_new[i] = ρ[i] - (dt/(2dr)) * Δρ

            Δm = (mom_flux[i+1] - mom_flux[i-1]) / r[i]^2
            v_new[i] = (ρ[i]*v[i] - (dt/(2dr))*Δm - ρ[i]*g_acc[i]) / ρ_new[i]
        end
    end

    # enforce Neumann BC at center
    ρ_new[1] = ρ_new[2]

    return ρ_new, v_new
end