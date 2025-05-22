




function get_initial_density(n, M, zones, interp_config)

    radius, dens = get_LE_density_profile(n, M)
    # Prepare for interpolation to lower-res grid
    dens_interp(var) = interpolation_dispatch(radius, dens, var, interp_config=interp_config)
    
    r = range(0.0, radius[end], length=zones)
    dr = r[2] - r[1]
    left_ghost_bound = r[1] - ghost_zones * dr
    total_zones = zones + 2 * ghost_zones
    rads = zeros(total_zones)
    for i in 1:total_zones
        rads[i] = left_ghost_bound + (i-1) * dr
    end

    # Interpolate densities directly
    rho = max.(dens_interp.(r), 0.0)  # Use broadcasting to apply condition
    
    mass_iters = length(r)
    mass = zeros(Float64, mass_iters)  # Pre-allocate mass array
    integrand = 4π * rho .* r.^2 #* dr
    dM = zeros(Float64, mass_iters)

    for i in 1:mass_iters-1
        interval = [r[i], r[i+1]]
        argument = [integrand[i], integrand[i+1]]
        dM[i] = gauss_quad(interval, argument, 4)
    end

    # Compute mass array via cumulative sum
    mass[2:end] .= cumsum(dM[2:end])

    return r, dr, rho, mass, dM, rads
end

function get_initial_velocity(mass, r)
    velocity = zeros(length(r))
    velocity[2:end] .= -sqrt.(2 .* G .* mass[2:end] ./ r[2:end])
    return velocity
end

function get_gravity_rhs(n, ρ, c_end; Φ_out=0.0)
    # 1) sample the source 4πGρ(r)
    rhs = 4π * G .* ρ 

    # 2) Neumann at r=0 via ghost → row 1 homogeneous
    rhs[1] = 0.0

    # 3) Dirichlet at outer radius: subtract c[n]*Φ_out
    rhs[n] -= c_end * Φ_out

    return rhs
end

function get_gravity(zones, dr, rho)
    # Build the spherical Poisson matrix A and the coefficient vector c
    A, c = build_spherical_poisson_matrix(zones, dr)
    
    # Compute the right-hand side (rhs) for the gravity equation
    rhs = get_gravity_rhs(zones, rho, c[end]; Φ_out=0.0)
    
    # Solve the system A * gravity = rhs using the backslash operator
    gravity = A \ rhs

    # Efficiently compute gravitational acceleration by vectorizing the diff
    iters = length(gravity)
    grav_accel = zeros(iters)
    for i in 2:iters-1
        grav_accel[i] = -(gravity[i+1] - gravity[i-1]) / (2 * dr)
    end
    grav_accel[1] = 0.0
    grav_accel[end] = -(gravity[end] - gravity[end-1]) / dr

    return gravity, grav_accel
end


function get_time_params(rho, velocity, dr)
    t_ff = sqrt(3π/(32 * G * maximum(rho)))
    dt = 0.5 * dr / maximum(abs.(velocity))
    return t_ff, dt
end


#=================================================================================================#
#                                                                                                 #
# To create an initial grid we need a minimized list of input variables. How can I create a grid  #
# with as little entries as possible?                                                             #
#                                                                                                 #
# Parameters:                                                                                     #
#                   i) zones         -> Number of grid points                                     #
#                  ii) M             -> Desired mass of polytrope                                 #
#                 iii) R             -> Desired radius of polytrope                               #
#                  iv) n             -> Adiabatic index (optional, default=5/3)                   #
#                   v) interp_config -> Interpolation Method (optional, defaults to cubic_spline) #
#                                                                                                 #
#=================================================================================================#

function create_initial_grid(n, M, zones, interp_config)
    r, dr, dens, mass, dM, rads = get_initial_density(n, M, zones, interp_config)
    velocity = get_initial_velocity(mass, r)
    gravity, grav_accel = get_gravity(zones, dr, dens)
    t_ff, dt = get_time_params(dens, velocity, dr)

    #=
    ghost_zones = 3  # You need to define this
    N = length(rads)
    radius = rads  # Assuming you meant this as the ghost radius grid

    ghost_rho = zeros(N)
    ghost_vel = zeros(N)
    ghost_grav_pot = zeros(N)
    ghost_grav_accel = zeros(N)

    # Apply symmetric boundary conditions
    ghost_rho[1:ghost_zones] .= reverse(dens[1:ghost_zones])
    ghost_rho[ghost_zones+1:end-ghost_zones] .= dens[ghost_zones+1:end-ghost_zones]
    ghost_rho[end-ghost_zones+1:end] .= reverse(dens[end-ghost_zones+1:end])

    ghost_vel[1:ghost_zones] .= reverse(velocity[1:ghost_zones])
    ghost_vel[ghost_zones+1:end-ghost_zones] .= velocity[ghost_zones+1:end-ghost_zones]
    ghost_vel[end-ghost_zones+1:end] .= reverse(velocity[end-ghost_zones+1:end])

    ghost_grav_pot[1:ghost_zones] .= reverse(gravity[1:ghost_zones])
    ghost_grav_pot[ghost_zones+1:end-ghost_zones] .= gravity[ghost_zones+1:end-ghost_zones]
    ghost_grav_pot[end-ghost_zones+1:end] .= reverse(gravity[end-ghost_zones+1:end])

    ghost_grav_accel[1:ghost_zones] .= reverse(grav_accel[1:ghost_zones])
    ghost_grav_accel[ghost_zones+1:end-ghost_zones] .= grav_accel[ghost_zones+1:end-ghost_zones]
    ghost_grav_accel[end-ghost_zones+1:end] .= reverse(grav_accel[end-ghost_zones+1:end])

    ghost_grid = [radius, ghost_rho, ghost_vel, ghost_grav_pot, ghost_grav_accel]
    =#
    return r, dr, dens, mass, dM, velocity, gravity, grav_accel, t_ff, dt
end
