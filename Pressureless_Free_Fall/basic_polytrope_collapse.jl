###############################################################################
# E.O.N.S. 1D Spherical Pressureless Collapse Solver (Physical Units)
# Refactored: no normalization, correct fluxes, modular functions
###############################################################################

#-- Dependencies ------------------------------------------------------------
include("../lane_emden.jl")
include("../../ARC/Numerical_Methods/interpolation/interpolate.jl")
include("../../ARC/Numerical_Methods/integration/methods/spatial/gaussian_quadrature.jl")
include("../Linear_Algebra/construct_tridiagonal_matrix.jl")
include("../Linear_Algebra/spherical_poisson_coefficients.jl")
include("build_initial_grid.jl")
using HDF5, Plots

#-- Physical Constants ------------------------------------------------------
const G     = 6.674e-8     # cm^3 / g / s^2
const M_sun = 1.98847e33   # g
const R_sun = 6.957e10     # cm

#-- Initial Grid (Physical Units) ------------------------------------------
zones       = 1000
ghost_zones = 6
n_polytrope = 5/3
M_star      = 10.0          # in solar masses

interp_cfg  = setup_interpolation(:cubic_spline, :optimize, false)
# Initial Quantities
r, dr, rho0, mass0, dM0, vel0, Φ0, g0, t_ff0, dt0 = create_initial_grid(n_polytrope, M_star, zones, interp_cfg)

iters = length(r)
rho = zeros(iters)
vel = zeros(iters)
Φ  = zeros(iters)
g  = zeros(iters)
dt = 0.5 * dr / dt0

function run_collapse_simulation!(
    rho::Vector{Float64},
    vel::Vector{Float64},
    g::Vector{Float64},
    dt::Float64,
    steps::Int,
    r::Vector{Float64},
    dr::Float64,
    zones::Int,
    G::Float64
)
    # Preallocate working arrays
    ρ_new = similar(rho)
    vel_new = similar(vel)

    # Build and factorize Poisson matrix once
    A, c = build_spherical_poisson_matrix(zones, dr)
    factor = lu(A)

    t = 0.0
    t_ff = sqrt(3π / (32G * maximum(rho)))
    time_history = Float64[]

    for step in 1:steps
        
        if t >= t_ff
            break
        end

        # Continuity and momentum updates
        @inbounds for i in 2:zones-1
            # continuity flux
            fluxρ = (rho[i+1]*vel[i+1]*r[i+1]^2 - rho[i-1]*vel[i-1]*r[i-1]^2) / (2dr)
            ρ = rho[i] - dt * fluxρ / r[i]^2
            ρ_new[i] = max(ρ, 0.0)

            # momentum flux
            fluxρv = (rho[i+1]*vel[i+1]^2*r[i+1]^2 - rho[i-1]*vel[i-1]^2*r[i-1]^2) / (2dr)
            ρv = rho[i]*vel[i] - dt * fluxρv / r[i]^2 + dt * rho[i] * g[i]

            vel_new[i] = ρ_new[i] == 0.0 ? 0.0 : ρv / ρ_new[i]
        end

        # BCs
        r[1] = 100000.0 # cm ~ 1 km
        ρ_new[1] = rho[1] - dt * (rho[2]*vel[2]*r[2]^2 - rho[1]*vel[1]*r[1]^2) / (dr*r[1]^2)
        r[1] = 0.0
        ρ_new[end] = 0.0
        vel_new[1] = 0.0
        vel_new[end] = 0.0

        # Update gravity: solve Poisson
        rhs = get_gravity_rhs(zones, ρ_new, c[end]; Φ_out=0.0)
        Φ = factor \ rhs

        g[1] = 0.0
        @inbounds for i in 2:zones-1
            g[i] = -(Φ[i+1] - Φ[i-1]) / (2dr)
        end
        g[end] = -(Φ[end] - Φ[end-1]) / dr

        # Advance time and adapt dt
        t += dt
        push!(time_history, t)
        t_ff = sqrt(3π / (32G * maximum(ρ_new)))
        dt = 0.5 * dr / maximum(abs.(vel_new))
        println(t)
        # Swap arrays
        rho .= ρ_new
        vel .= vel_new
    end

    return rho, vel, g, time_history
end


steps = 100
rho0, vel0, g0, time = run_collapse_simulation!(rho0, vel0, g0, dt0, steps, collect(r), dr, zones, G)

println(sum(time))
#plot(r, rho0, label="Density", xlabel="Radius (cm)", ylabel="Density (g/cm^3)", title="Density Profile", legend=:topright)
r = collect(r)
r[1] = 1e-10
rho0_clipped = @. max(rho0, 1e-1)

plot(r, rho0_clipped,
    label = "Density",
    xlabel = "Radius (cm)",
    ylabel = "Density (g/cm³)",
    title = "Density Profile",
    legend = :topright,
    yscale = :log10,
    minorgrid = true,
    minorgridalpha = 0.1,
    grid = true
)

#=
#plot(r, vel0, label="Velocity", xlabel="Radius (cm)", ylabel="Velocity (cm/s)", title="Velocity Profile", legend=:topright)
plot(r, g0,
    label = "Gravitational Acceleration",
    xlabel = "Radius (cm)",
    ylabel = "Gravitational Acceleration (cm/s²)",
    title = "Gravitational Acceleration Profile",
    legend = :topright,
    grid = true
)
=#
