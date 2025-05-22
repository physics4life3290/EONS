




function evolve_collapse(rho0, velocity0, gravity0, grav_accel0, r, dr, t_ff, dt;
    steps=300, save_every=10, output_file="collapse_data.h5")


    zones = length(r)
    rho = copy(rho0)
    velocity = copy(velocity0)
    gravity = copy(gravity0)
    grav_accel = copy(grav_accel0)
    ff_timescale = [t_ff, ]

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

        rho, velocity = update_free_fall!(rho, velocity, r, grav_accel, dr, dt)
        #rho .= update_density!(rho, velocity, r, dt, dr)
    
        A, c = build_spherical_poisson_matrix(zones, dr)
        rhs = get_gravity_rhs(zones, rho, c[end]; Φ_out=0.0)
        
        gravity = A \ rhs
        grav_accel = vcat(0.0, -diff(gravity))

        #velocity .= update_velocity!(velocity, grav_accel, dt, dr)
        
        t_ff = sqrt(3π/(32 * G * maximum(rho)))
        push!(ff_timescale, t_ff) 
        
        # save every `save_every` steps
        if step % save_every == 0
            groupname = "/snapshots/step_$(step)"
            g = create_group(h5file, groupname)

            g["radius"] = collect(r)
            g["density"] = rho
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
    
    #return  collect(r), rho, velocity, gravity, grav_accel, ff_timescale, dt
end
