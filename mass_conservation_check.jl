




function mass_conservation_check(time, rho, velocity, rho_new, velocity_new)
    time_interval = [time[1], time[2]] 
    outer_interval = [rho[end]*velocity[end], rho_new[end]*velocity_new[end]]
    inner_interval = [rho[1]*velocity[1], rho_new[1]*velocity_new[1]]
    gauss_quad(time_interval, outer_interval, 5) - gauss_quad(time_interval, inner_interval, 5)


    # Compute initial mass
    integrand_initial = rho .* (4π .* r.^2)
    mass_initial = gauss_quad(r, integrand_initial, 4)

    # Compute final mass
    integrand_final = rho_new .* (4π .* r_new.^2)
    mass_final = gauss_quad(r_new, integrand_final, 4)

    # Difference (should be close to 0)
    mass_difference = mass_final - mass_initial
    mass_relative_error = abs(mass_difference) / mass_initial

    println("Initial mass: $mass_initial")
    println("Final mass: $mass_final")
    println("Relative error: $mass_relative_error")
end
