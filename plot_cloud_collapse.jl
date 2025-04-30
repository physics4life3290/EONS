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
                        every::Int=1, title_prefix="")
    steps = sort(collect(keys(snapshots)))
    plt = plot(title="$title_prefix $field", xlabel="Radius (cm)", ylabel=field)
    for (i, step) in enumerate(steps)
        if i % every == 0
            snap = snapshots[step]
            plot!(plt, snap["radius"], snap[field], label="t=$(round(snap["time"], sigdigits=3))s")
        end
    end
    return plt
end

# Usage example
h5_filename = "collapse_run.h5"
snapshots = load_snapshots(h5_filename)

# Plot snapshots for different quantities
p1 = plot_snapshots(snapshots, "density", every=3, title_prefix="Collapse:")
p2 = plot_snapshots(snapshots, "velocity", every=3, title_prefix="Collapse:")
p3 = plot_snapshots(snapshots, "gravity", every=3, title_prefix="Collapse:")
p4 = plot_snapshots(snapshots, "grav_accel", every=3, title_prefix="Collapse:")

# Display the plots
#plot(p1)
#plot(p2)
#plot(p3)
plot(p4)
