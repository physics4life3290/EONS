using HDF5, Plots

function load_snapshots(h5_filename::String)
    snapshots = Dict{Int, Dict}()

    h5file = h5open(h5_filename, "r")
    group = h5file["/snapshots"]

    for name in keys(group)
        step = parse(Int, split(name, "_")[end])
        snapshot = Dict(
            "radius" => read(group[name]["radius"]),
            "density" => read(group[name]["density"]),
            "pressure" => read(group[name]["pressure"]),
            #"temperature" => read(group[name]["temperature"]),
            #"energy" => read(group[name]["energy"]),
            "velocity" => read(group[name]["velocity"]),
            "gravity" => read(group[name]["gravity"]),
            "grav_accel" => read(group[name]["grav_accel"]),
            "time" => read(group[name]["time"]),
            "dt" => read(group[name]["dt"]),
        )
        snapshots[step] = snapshot
    end

    close(h5file)
    return snapshots
end

function plot_snapshots(snapshots::Dict{Int, Dict}, field::String;
                        every::Int=1, title_prefix="",
                        log_scale_x::Bool=false, log_scale_y::Bool=false,
                        xlims::Tuple{Float64, Float64}=(nothing, nothing))

    steps = sort(collect(keys(snapshots)))
    plt = plot(title="$title_prefix $field", xlabel="Radius (cm)", ylabel=field)
    
    for (i, step) in enumerate(steps)
        if i % every == 0
            snap = snapshots[step]
            x_data = snap["radius"]
            y_data = snap[field]
            
            # Apply log scale to the x-axis (radius)
            if log_scale_x
                x_data = log10.(x_data .+ 1e-10)  # Adding a small value to avoid log(0)
            end
            
            # Apply log scale to the y-axis (density, velocity, etc.)
            if log_scale_y && field == "density"
                y_data = log10.(y_data .+ 1e-10)  # Adding a small value to avoid log(0)
            end

            # Apply log scale to the y-axis (density, velocity, etc.)
            if log_scale_y && field == "pressure" 
                y_data = log10.(y_data .+ 1e-10)  # Adding a small value to avoid log(0)
            end

            plot!(plt, x_data, y_data, label="t=$(round(snap["time"], sigdigits=3))s")
        end
    end
    
    # Set xlims if specified
    if xlims != (nothing, nothing)
        xlims!(plt, xlims)
    end


    return plt
end

R_sun = 6.957E10

# Usage example
h5_filename = "isothermal_collapse/isothermal_collapse_run.h5"
snapshots = load_snapshots(h5_filename)

# Plot snapshots for different quantities with custom xlims
p1 = plot_snapshots(snapshots, "density", every=1, title_prefix="Collapse:", log_scale_x=false, log_scale_y=false, xlims=(0.0, R_sun))
p2 = plot_snapshots(snapshots, "velocity", every=1, title_prefix="Collapse:", xlims=(0.0, R_sun))
p3 = plot_snapshots(snapshots, "gravity", every=1, title_prefix="Collapse:", xlims=(0.0, R_sun))
p4 = plot_snapshots(snapshots, "grav_accel", every=1, title_prefix="Collapse:", xlims=(0.0, R_sun))
p5 = plot_snapshots(snapshots, "pressure", every=1, title_prefix="Collapse:", xlims=(0.0, R_sun))
#p6 = plot_snapshots(snapshots, "temperature", every=1, title_prefix="Collapse:", xlims=(0.0, R_sun))
#p7 = plot_snapshots(snapshots, "energy", every=1, title_prefix="Collapse:", log_scale_y = true , xlims=(0.0, R_sun))

# Display the plots
#plot(p1)
plot(p2)
#plot(p3)
#plot(p4)
#plot(p5)
#plot(p6)
#plot(p7)


