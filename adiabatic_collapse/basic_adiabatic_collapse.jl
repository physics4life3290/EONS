











include("../lane_emden.jl")
include("../../ARC/Numerical_Methods/interpolation/interpolate.jl")
include("../../ARC/Numerical_Methods/integration/methods/spatial/gaussian_quadrature.jl")
include("../Linear_Algebra/construct_tridiagonal_matrix.jl")
include("../Linear_Algebra/spherical_poisson_coefficients.jl")
include("build_initial_grid.jl")
include("evolve_adiabatic.jl")
include("adiabatic_finite_update.jl")

using HDF5, Plots

n = 5/3
M = 1.0
zones = 1000

G = 6.674E-8 # Gravitational constant in cm^3/g/s^2
k_B = 1.380649E-16  # erg ⋅ K^-1
mu = 1.00784
m_p = 1.6726219E-24 # g
R_sun = 6.957E10 # solar radius in cm
M_sun = 1.98847E33 # g
rad_cnst = 7.5646E-15 # erg cm^-3 K^-4
R_gas_cnst = 8.314462618E7 # erg mol^-1 K^-1
time = [0.0,]

method = :cubic_spline 
mode = :optimize 
extrapolate = false
interp_config = setup_interpolation(method, mode, extrapolate)

r, dr, rho, mass, dM, pressure, temperature, velocity, gravity, grav_accel, energy, t_ff, dt = create_initial_grid(n, M, zones, interp_config)

evolve_collapse(rho ,velocity, gravity, grav_accel, pressure, energy, r, dr, t_ff, dt, steps=200, save_every=10, output_file="adiabatic_collapse_run.h5")
