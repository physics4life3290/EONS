




function get_initial_density(n, M, zones, interp_config)

    radius, dens = get_LE_density_profile(n, M)
    # Prepare for interpolation to lower-res grid
    dens_interp(var) = interpolation_dispatch(radius, dens, var, interp_config=interp_config)
    
    r = range(0, radius[end], length=zones)
    dr = r[2] - r[1]

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

    return r, dr, rho, mass, dM
end

function get_initial_velocity(mass, r)
    velocity = -sqrt.(2 * G .* mass ./ r)
    velocity[1] = 0.0
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
    grav_accel = -diff(gravity)  # Use `diff` to compute the derivative
    grav_accel = vcat(0.0, grav_accel)  # Add zero for the first element to match dimensions

    return gravity, grav_accel
end


function get_time_params(rho, velocity, dr)
    t_ff = sqrt(3π/(32 * G * maximum(rho)))
    dt = 0.5 * dr / maximum(abs.(velocity))
    #dt = dr / 3E10
    return t_ff, dt
end

function create_initial_grid(n, M, zones, interp_config)

    r, dr, rho, mass, dM = get_initial_density(n, M, zones, interp_config)
    #velocity = zeros(length(r))
    velocity = get_initial_velocity(dM, r)
    gravity, grav_accel = get_gravity(zones, dr, rho)
    
    total_mass = sum(dM)
    P_c = (3 * G / 8π) * (total_mass^2 / r[end]^4)
    T_c = (P_c * mu * m_p) / (rho[1] * k_B)

    pressure = zeros(zones)
    temperature = zeros(zones)
    pressure[1] = P_c
    temperature[1] = T_c

    for i in 2:zones
        pressure[i] = pressure[i-1] * (rho[i]/rho[i-1])^(n)
        temperature[i] = temperature[i-1] * (rho[i]/rho[i-1])^(n-1)
    end
    
    energy = 0.5 .* rho .* velocity.^2 .+ rho .* (R_gas_cnst / (n - 1)) .* temperature

    t_ff, dt = get_time_params(rho, velocity, dr)

    return r, dr, rho, mass, dM, pressure, temperature, velocity, gravity, grav_accel, energy, t_ff, dt
end