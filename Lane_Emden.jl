function lane_emden(n::Float64, ξ_max::Float64, h::Float64)
    # Initial conditions near ξ = 0
    θ₀ = 1.0
    y₀ = 0.0
    ξ₀ = 1e-6

    ξ_vals = [0.0, ξ₀]
    θ_vals = [θ₀, θ₀]
    y_vals = [y₀, y₀]

    ξ = ξ₀
    θ = θ₀
    y = y₀

    # Define RHS of the system: dθ/dξ = y, dy/dξ = -θ^n - (2/ξ)y
    dθ(ξ, θ, y) = y
    dy(ξ, θ, y) = (θ > 0 ? -θ^n : 0.0) - (2 / ξ) * y

    while ξ < ξ_max
        if θ <= 0
            @warn "θ has become non-positive, stopping integration."
            break
        end

        # RK4 integration
        k1_θ = h * dθ(ξ, θ, y)
        k1_y = h * dy(ξ, θ, y)

        k2_θ = h * dθ(ξ + h/2, θ + k1_θ/2, y + k1_y/2)
        k2_y = h * dy(ξ + h/2, θ + k1_θ/2, y + k1_y/2)

        k3_θ = h * dθ(ξ + h/2, θ + k2_θ/2, y + k2_y/2)
        k3_y = h * dy(ξ + h/2, θ + k2_θ/2, y + k2_y/2)

        k4_θ = h * dθ(ξ + h, θ + k3_θ, y + k3_y)
        k4_y = h * dy(ξ + h, θ + k3_θ, y + k3_y)

        θ_new = θ + (k1_θ + 2*k2_θ + 2*k3_θ + k4_θ)/6
        y_new = y + (k1_y + 2*k2_y + 2*k3_y + k4_y)/6
        ξ += h

        if θ_new < 0
            break
        end

        push!(ξ_vals, ξ)
        push!(θ_vals, θ_new)
        push!(y_vals, y_new)

        θ = θ_new
        y = y_new
    end

    return ξ_vals, θ_vals, y_vals
end




#=
α = 1 # Scale factor for radius
ρ_c = 5.7
n = 1.0#5/3
R = 6.957E10 # 1 solar radius in cm

# Example usage:
ξ_vals, θ_vals, y_vals = lane_emden(n, 10.0, 1E-3)

ξ_vals = ξ_vals ./ maximum(ξ_vals)
# Output radius where θ crosses zero
iters = length(θ_vals)
for i in 2:iters
    if θ_vals[i] < 0
        println("Approximate surface radius (ξ where θ ≈ 0): ", ξ_vals[i-1])
        break
    end
end

radius = α .* ξ_vals .* R 
density = ρ_c .* θ_vals .^ n

counter = 0
iters = length(radius)
for i in 1:iters
    if density[i] > ρ_c
        density[i] = ρ_c
        global counter += 1
    end
end

dρ = zeros(iters)
dρ[1] = 0.0
dρ[2:end] = diff(density)

density_stiffness = dρ ./ density

using Plots
#plot(ξ_vals, θ_vals, label="Theta Profile")
plot(radius[3:end], abs.(density_stiffness[3:end]), label="Density Profile", yscale=:log10)
#plot(radius, 1e-3 ./ radius, yscale=:log10)
#plot(radius, density, label="Density Profile")
#println(length(density))
#println(radius[end])
#println(density[end])
#println(density[2])
#println(radius[1:10])
=#