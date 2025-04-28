



R_sun = 6.957E10 # solar radius in cm
R = R_sun
G = 6.674E-8 # Gravitational constant in cm^3/g/s^2
zones = 1000
rho_c = 5.7
r = collect(range(0, R, zones))
dr = r[2] - r[1]
rho = rho_c .* (1 .- (r ./ R) .^ 2)

mass_iters = length(r)
mass = zeros(mass_iters)
dM = zeros(mass_iters)

for i in 1:mass_iters-1
    dM[i+1] = 4π * rho[i] * r[i]^2 * dr
    mass[i+1] = mass[i] + dM[i+1]
end

velocity = -sqrt.(2 * G .* mass ./ r)
velocity[1] = 0.0

function build_tridiagonal(a::Vector, b::Vector, c::Vector)
    n = length(b)
    T = zeros(Float64, n, n)
    for i in 1:n
        T[i,i] = b[i]
        if i > 1
            T[i,i-1] = a[i]
        end
        if i < n
            T[i,i+1] = c[i]
        end
    end
    return T
end

function build_spherical_tridiagonal_matrix(n, Δr)
    a = zeros(n); b = zeros(n); c = zeros(n)
    for i in 1:n
        if i == 1
            # enforce f'(0)=0 via ghost: 
            b[i] = -2/Δr^2
            c[i] = +2/Δr^2
        else
            r = i*Δr
            a[i] = 1/Δr^2 - 1/(r*Δr)   # sub-diagonal
            b[i] = -2/Δr^2             # diagonal
            c[i] = 1/Δr^2 + 1/(r*Δr)   # super-diagonal
        end
    end
    return build_tridiagonal(a,b,c), c
end

function make_rhs(n, ρ, G, c_end; Φ_out=0.0)
    # 1) sample the source 4πGρ(r)
    rhs = 4π * G .* ρ 

    # 2) Neumann at r=0 via ghost → row 1 homogeneous
    rhs[1] = 0.0

    # 3) Dirichlet at outer radius: subtract c[n]*Φ_out
    rhs[n] -= c_end * Φ_out

    return rhs
end

# build a, b, c as before (flux form or central diffs)
A, c = build_spherical_tridiagonal_matrix(zones, dr)
rhs = make_rhs(zones, rho, G, c[end]; Φ_out=0.0)

gravity = A \ rhs

grav_accel = zeros(length(rhs))
grav_accel[2:end] = -diff(gravity)

t_ff = sqrt(3π/(32 * G * maximum(rho)))
println(t_ff)
dt = 0.5 * dr / maximum(abs.(velocity))


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
        ρ_new[i-1] = ρ[i] - (Δt/(2Δr)) * (1/r[i]^2) * (f[i+1] - f[i-1])
    end

    return ρ_new
end

function update_velocity!(
    v::AbstractVector,
    g::AbstractVector,
    Δt::Float64,
    Δr::Float64
)
    N = length(v)
    v_new = zeros(N)

    @inbounds for i in 2:N-1
        v_new[i-1] = v[i-1] + g[i-1]*Δt - (v[i-1]*Δt)/Δr * (v[i+1] - v[i-1])
    end

    return v_new
end


time = [0.0,]

for step in 1:250
    # recompute CFL dt, but guard zero-velocity
    push!(time, dt)
    total_time = sum(time)
    
    vmax = maximum(abs.(velocity))
    if vmax > 0
      global dt = 0.5 * dr /  vmax
    end
  
  
    # update density and velocity
    global rho      .= update_density!(rho, velocity, r, dt, dr)
    global velocity .= update_velocity!(velocity, grav_accel, dt, dr)
    
    iters = length(rho)
    for i in 1:iters
        if rho[i] < 0.0
            rho[i] = 0.0
        end
    end
    # recompute gravity
    global A, c    = build_spherical_tridiagonal_matrix(zones, dr)
    global rhs     = make_rhs(zones, rho, G, c[end]; Φ_out=0.0)
    global gravity = A \ rhs
    global grav_accel = vcat(0.0, -diff(gravity))
  
    # (optionally print progress)
    if step % 100 == 0
        println(" step=", step,
              "  maxρ=", round(maximum(rho), digits=4),
              "  dt=", round(dt, sigdigits=3))
        println("The total time passed in the simulation is: $(total_time)")
    end
  end


using Plots





plot(r, velocity ./ maximum(abs.(velocity)))
savefig("velocity_plot.pdf")


plot(r, gravity, xlabel="r", ylabel="Gravitational Potential", label="Gravitational Potential")
savefig("gravity_plot.pdf")

plot(r, grav_accel)
savefig("grav_accel.pdf")

plot(r, rho, label="Evolved Density")
savefig("evolved_density.pdf")

plot(r, mass ./ maximum(mass))
plot!(r, dM ./ maximum(dM))
savefig("mass_plot.pdf")
