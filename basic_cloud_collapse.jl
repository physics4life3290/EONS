




include("Lane_Emden.jl")
include("Cubic_spline.jl")
include("Linear_Algebra/construct_tridiagonal_matrix.jl")
include("Linear_Algebra/spherical_poisson_coefficients.jl")
include("build_initial_grid.jl")
include("free_fall_finite_diff_update.jl")
include("evolve_free_fall.jl")
using HDF5, Plots

R_sun = 6.957E10 # solar radius in cm
R = R_sun
G = 6.674E-8 # Gravitational constant in cm^3/g/s^2
zones = 1000
rho_c = 5.7
α=1.0
n = 5/3
time = [0.0,]

r, dr, rho, mass, dM, velocity, gravity, grav_accel, t_ff, dt = create_initial_grid(n, zones, R, rho_c)

evolve_collapse(rho, velocity, gravity, grav_accel, r, dr, t_ff, dt, steps=280, save_every=1, output_file="collapse_run.h5")