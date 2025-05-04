




function cubic_spline(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x)
    α = [0.0; 3.0 * (y[3:end] .- y[2:end-1]) ./ h[2:end] .- 3.0 * (y[2:end-1] .- y[1:end-2]) ./ h[1:end-1]]

    # Tridiagonal system
    l = ones(n)
    μ = zeros(n)
    z = zeros(n)

    for i in 2:n-1
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * μ[i-1]
        μ[i] = h[i] / l[i]
        z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
    end

    # Back substitution
    b = zeros(n-1)
    c = zeros(n)
    d = zeros(n-1)

    for j in n-1:-1:1
        c[j] = z[j] - μ[j] * c[j+1]
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2*c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3*h[j])
    end

    # Create interpolant function
    function spline_interp(xi)
        i = findfirst(xi .< x[2:end])  # Find appropriate interval
        i = isnothing(i) ? n-1 : i     # Clamp to last interval if out of bounds

        dx = xi - x[i]
        return y[i] + b[i]*dx + c[i]*dx^2 + d[i]*dx^3
    end

    return spline_interp
end
#=
x = collect(0.0:0.1:10.0)
y = sin.(x)
spline = cubic_spline(x, y)

# Single interpolation
println(spline(1.75))

x = collect(range(0.0, 100, 1001))
f(var) = var^2
y = f.(x)
spline = cubic_spline(x, y)
println(spline(1.75))  # interpolated value at x = 2.5
=#
