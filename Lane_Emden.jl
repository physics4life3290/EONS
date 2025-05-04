

# Function for stellar radius
function stellar_radius(M::Float64; units::Bool=true)
    if M < 0.1
        error("Mass too low for main-sequence approximation (< 0.1 M_sun).")
    elseif M <= 1.5
        R = M^0.8
    elseif M <= 10
        R = M^0.57
    else
        R = M^0.35
    end
    return units ? R : R * 6.957e10  # Convert to cm if needed (1 R_sun = 6.957e10 cm)
end

# Lane-Emden solver
function solve_lane_emden(n, xi_max; h=0.001)
    ε = h
    ξ = [0.0, ε]
    θ = [1.0, 1.0 - ε^2/6]
    ϕ = [0.0, -ε/3]

    function f(x, θ, ϕ)
        θ_eff = max(θ, 0.0)
        dθ = ϕ
        dϕ = -2/x*ϕ - θ_eff^n
        return dθ, dϕ
    end

    while ξ[end] < xi_max && θ[end] > 0
        x, y, y′ = ξ[end], θ[end], ϕ[end]
        k1y,  k1y′  = f(x,           y,            y′)
        k2y,  k2y′  = f(x + h/2,     y + k1y*h/2,   y′ + k1y′*h/2)
        k3y,  k3y′  = f(x + h/2,     y + k2y*h/2,   y′ + k2y′*h/2)
        k4y,  k4y′  = f(x + h,       y + k3y*h,     y′ + k3y′*h)

        θ_next = y  + (h/6)*(k1y  + 2k2y  + 2k3y  + k4y)
        ϕ_next = y′ + (h/6)*(k1y′ + 2k2y′ + 2k3y′ + k4y′)
        x_next = x + h

        push!(ξ, x_next)
        push!(θ, θ_next)
        push!(ϕ, ϕ_next)
    end

    return ξ, θ, ϕ
end

# Function to find the first zero
function first_zero(ξ, θ)
    for i in 1:length(θ)-1
        if θ[i] > 0 && θ[i+1] <= 0
            return ξ[i] - θ[i]*(ξ[i+1]-ξ[i])/(θ[i+1]-θ[i])
        end
    end
    return NaN
end

# Function for generating density profile
function get_LE_density_profile(n, M)
    ξ, θ, ϕ = solve_lane_emden(n, 10.0)
    xi1 = first_zero(ξ, θ)
    R = stellar_radius(M, units=false)

    alpha = R / xi1

    dtheta_dxi_surface = ϕ[findfirst(t -> t <= 0, θ) - 1]
    mass_prefactor = -4π * alpha^3 * xi1^2 * dtheta_dxi_surface
    rho_c = M * M_sun / mass_prefactor

    N = findfirst(t -> t <= 0, θ)
    if N !== nothing
        ξ = ξ[1:N]
        θ = θ[1:N]
        r = alpha .* ξ
    else
        r = alpha .* ξ
    end

    ρ = rho_c .* (max.(θ,0.0) .^ n)

    return r, ρ
end

