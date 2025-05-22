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
zones       = 100
ghost_zones = 6
n_polytrope = 5/3
M_star      = 1.0          # in solar masses

interp_cfg  = setup_interpolation(:cubic_spline, :optimize, false)
# Initial Quantities
r, dr, rho0, mass0, dM0, vel0, Φ0, g0, t_ff0, dt0 = create_initial_grid(n_polytrope, M_star, zones, interp_cfg)

iters = length(r)
rho = zeros(iters)
vel = zeros(iters)
Φ  = zeros(iters)
g  = zeros(iters)
dt = 0.5 * dr / 3e10 #dt0
t_ff = t_ff0

function run_collapse_simulation!(rho0, vel0, g0, dt0, steps, r, dr, zones, G)
    rho = copy(rho0)
    vel = copy(vel0)
    time = [0.0,]

    for time_step in 1:steps
        for i in 2:zones-1
            # Continuity equation
            flux_rho = (rho0[i+1]*vel0[i+1]*r[i+1]^2 - rho0[i-1]*vel0[i-1]*r[i-1]^2) / (2dr)
            rho[i] = rho0[i] - dt0 * flux_rho / r[i]^2
            #rho[i] = ((rho0[i+1]+rho0[i-1])/2) - dt0 * flux_rho / r[i]^2
            rho[i] = max(rho[i], 0.0)

            # Momentum equation (update rho*v first)
            flux_rho_v = (rho0[i+1]*vel0[i+1]^2*r[i+1]^2 - rho0[i-1]*vel0[i-1]^2*r[i-1]^2) / (2dr)
            rho_v = rho0[i]*vel0[i] - dt0 * flux_rho_v / r[i]^2 + dt0 * rho0[i] * g0[i]
            #rho_v = (rho0[i+1]*vel0[i+1] + rho0[i-1]*vel0[i-1])/2 - dt0 * flux_rho_v / r[i]^2 + dt0 * rho0[i] * g0[i]
            if rho[i] == 0.0
                vel[i] = 0.0
            else
                # Update velocity
                vel[i] = rho_v / rho[i]
            end
        end


        rho[1] = rho0[1] - dt0/(dr) * (rho0[2]*vel0[2]*r[2]^2 - rho0[1]*vel0[1]*r[1]^2) / (1000)^2
        rho[end] = 0.0
        vel[1] = vel[end] = 0.0

        Coeff_mat, c = build_spherical_poisson_matrix(zones, dr)
        rhs = get_gravity_rhs(zones, rho, c[end]; Φ_out=0.0)
        Φ = Coeff_mat \ rhs

        iters = length(r)
        g = zeros(iters)
        for i in 2:iters-1
            g[i] = -(Φ[i+1] - Φ[i-1]) / (2 * dr)
        end
        g[1] = 0.0
        g[end] = -(Φ[end] - Φ[end-1]) / dr

        dt = 0.5 * dr / maximum(abs.(vel))
        t_ff = sqrt(3π / (32 * G * maximum(rho)))

        push!(time, dt)
        rho0 .= rho
        vel0 .= vel
        g0 .= g
        dt0 = 0.5 * dr/ 3E10 #dt

    end

    return rho0, vel0, g0, time
end

steps = 100
rho0, vel0, g0, time = run_collapse_simulation!(rho0, vel0, g0, dt0, steps, r, dr, zones, G)

println(sum(time))
#plot(r, rho0, label="Density", xlabel="Radius (cm)", ylabel="Density (g/cm^3)", title="Density Profile", legend=:topright)
r = collect(r)
r[1] = 1e-10
rho0_clipped = @. max(rho0, 1e-10)

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