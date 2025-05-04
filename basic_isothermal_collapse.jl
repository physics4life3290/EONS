











include("lane_emden.jl")
include("../ARC/Numerical_Methods/interpolation/interpolate.jl")
include("../ARC/Numerical_Methods/integration/methods/spatial/gaussian_quadrature.jl")
include("Linear_Algebra/construct_tridiagonal_matrix.jl")
include("Linear_Algebra/spherical_poisson_coefficients.jl")
using HDF5, Plots

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
    return t_ff, dt
end

function create_initial_grid(n, M, zones, interp_config)

    r, dr, rho, mass, dM = get_initial_density(n, M, zones, interp_config)
    velocity = get_initial_velocity(mass, r)
    gravity, grav_accel = get_gravity(zones, dr, rho)
    
    pressure = rho .* k_B .* T /(mu .* m_p)
    t_ff, dt = get_time_params(rho, velocity, dr)

    return r, dr, rho, mass, dM, pressure, velocity, gravity, grav_accel, t_ff, dt
end

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

function evolve_collapse(rho0, velocity0, gravity0, grav_accel0, pressure0, r, dr, t_ff, dt;
                         steps=300, save_every=10, output_file="collapse_data.h5")

    zones = length(r)
    rho = copy(rho0)
    velocity = copy(velocity0)
    gravity = copy(gravity0)
    grav_accel = copy(grav_accel0)
    P = copy(pressure0)

    h5file = h5open(output_file, "w")

    # store some metadata
    h5file["/metadata/t_ff"] = t_ff
    h5file["/metadata/initial_dt"] = dt

    step = 0.0
    for step in 1:steps
    #while sum(time) < t_ff
        #step += 1.0
        push!(time, dt)
        total_time = sum(time)

        vmax = maximum(abs.(velocity))
        if vmax > 0
            dt = 0.5 * dr / vmax
        end

        rho .= update_density!(rho, velocity, r, dt, dr)
        rho_iters = length(rho)
        for i in 1:rho_iters
            if rho[i] < 0.0
                rho[i] = 0.0
            end
        end

        P = rho .* k_B  .* T ./ (mu .* m_p)

        A, c = build_spherical_poisson_matrix(zones, dr)
        rhs = get_gravity_rhs(zones, rho, c[end]; Φ_out=0.0)
        
        gravity = A \ rhs
        grav_accel = vcat(0.0, -diff(gravity))
        
        velocity .= update_velocity!(velocity, grav_accel, rho, P, dt, dr)
    
    
    
        # save every `save_every` steps
        if step % save_every == 0
            groupname = "/snapshots/step_$(step)"
            g = create_group(h5file, groupname)

            g["radius"] = collect(r)
            g["density"] = rho
            g["pressure"] = P
            g["velocity"] = velocity
            g["gravity"] = gravity
            g["grav_accel"] = grav_accel
            g["time"] = total_time
            g["dt"] = dt

            close(g)
        end
        
        # (optional) print progress
        if step % 100 == 0
            println("Step $step: max(ρ) = ", round(maximum(rho), digits=5), 
                    ", dt = ", round(dt, sigdigits=3),
                    ", total time = ", round(total_time, sigdigits=4))
        end
        
    end
    close(h5file)
    
    println("Simulation complete. Data saved to: $output_file")
end

n = 5/3
M = 1.0
zones = 1000
T = 30 # Kelvin

G = 6.674E-8 # Gravitational constant in cm^3/g/s^2
k_B = 1.380649E-16  # erg ⋅ K^-1
mu = 1.00784
m_p = 1.6726219E-24 # g
R_sun = 6.957E10 # solar radius in cm
M_sun = 1.98847E33

time = [0.0,]

method = :cubic_spline 
mode = :optimize 
extrapolate = false
interp_config = setup_interpolation(method, mode, extrapolate)

r, dr, rho, mass, dM, pressure, velocity, gravity, grav_accel, t_ff, dt = create_initial_grid(n, M, zones, interp_config)

evolve_collapse(rho ,velocity, gravity, grav_accel, pressure, r, dr, t_ff, dt, steps=700, save_every=10, output_file="isothermal_collapse_run.h5")


#=
#plot(r, rho)
plot(r, new_vel)
dP = zeros(length(pressure))
dP[2:end] = diff(pressure)

iters = length(rho)
pressure_accel = zeros(iters)
for i in 1:iters
    if rho[i] == 0.0
        pressure_accel[i] = 0.0
    else 
        pressure_accel[i] = dP[i] ./ rho[i]
    end
end

dP_dr = dP ./ dr
println(dP[4])

#plot(r, pressure ./ maximum(pressure), label="Pressure")
#plot!(r, dP_dr ./ maximum(abs.(dP_dr)), label="Pressure Gradient")
#plot(r, pressure_accel , label="Pressure Acceleration")
#plot!(r, gravity ./ maximum(abs.(gravity)), label="Gravity")
#plot!(r, -grav_accel, label="Gravitational Acceleration")

#plot(r[2:end-1], abs.(pressure_accel[2:end-1]), label="Pressure Gradient / Density", yscale=:log10)
#plot!(r, grav_accel, label="Gravitational Acceleration")
=#