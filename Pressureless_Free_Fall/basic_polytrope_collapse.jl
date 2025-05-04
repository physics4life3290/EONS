




include("../lane_emden.jl")
include("../../ARC/Numerical_Methods/interpolation/interpolate.jl")
include("../../ARC/Numerical_Methods/integration/methods/spatial/gaussian_quadrature.jl")
include("../Linear_Algebra/construct_tridiagonal_matrix.jl")
include("../Linear_Algebra/spherical_poisson_coefficients.jl")
include("build_initial_grid.jl")
include("free_fall_finite_diff_update.jl")
include("evolve_free_fall.jl")
using HDF5, Plots


G = 6.674E-8 # Gravitational constant in cm^3/g/s^2
M_sun = 1.98847E33
R_sun = 6.957E10

zones = 1000
n = 5/3
M = 1.0 # 1 Solar mass 

time = [0.0,]

method = :cubic_spline 
mode = :optimize 
extrapolate = false
interp_config = setup_interpolation(method, mode, extrapolate)

r, dr, rho, mass, dM, velocity, gravity, grav_accel, t_ff, dt = create_initial_grid(n, M, zones, interp_config)

evolve_collapse(rho, velocity, gravity, grav_accel, r, dr, t_ff, dt, steps=100, save_every=1, output_file="pressureless_collapse_run.h5")