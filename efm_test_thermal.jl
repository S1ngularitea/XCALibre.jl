using XCALibre
using CSV
# using CUDA # Uncomment to run on NVIDIA GPUs
# using AMDGPU # Uncomment to run on AMD GPUs

grids_dir = pkgdir(XCALibre, "Test_Meshes/");
grid = "initial_efm_mesh_CD1_2.unv";
mesh_file = joinpath(grids_dir, grid);

mesh = UNV2D_mesh(mesh_file, scale=0.001);

# Select backend and setup hardware
backend = CPU();
# backend = CUDABackend() # run on NVIDIA GPUs
# backend = ROCBackend() # run on AMD GPUs

hardware = Hardware(backend=backend, workgroup=1024);
# hardware = Hardware(backend=backend, workgroup=32) # use for GPU backends

mesh_dev = mesh; # use this line to run on CPU
# mesh_dev = adapt(backend, mesh)  # Uncomment to run on GPU 

begin # Model Info
inlet_width = 0.510; # m
inlet_height = 0.05; # m
inlet_area = inlet_width*inlet_height; # m^2

if isfile("Model_Input.csv")
    test_case = 7;
    input_parameters = CSV.File("Model_Input.csv"); #File containing the different test cases from paper "Modeling of Partially Wetting Liquid Film Using an Enhanced Thin Film Model for Aero-Engine Bearing Chamber Applications" by Kuldeep Singh et. al

    inlet_flow_rate = input_parameters.Q[test_case]; # m^3/s
    β = input_parameters.Beta[test_case]; # empisical value from paper
    σ = input_parameters.Sigma[test_case]
    θm = input_parameters.Theta_m[test_case]
    ϕ = input_parameters.Phi[test_case]
else
    inlet_flow_rate = 25.53e-5
    β = 3
    σ = 0.069
    θm = 75
    ϕ = 90
end

# EFM_DT sets the requested initial time step.
EFM_DT = parse(Float64, get(ENV, "EFM_DT", "5e-4"))
# EFM_ITERS sets the number of solver iterations.
EFM_ITERS = parse(Int, get(ENV, "EFM_ITERS", "20000"))
# EFM_WRITE_INTERVAL sets output frequency; use -1 to disable writes.
EFM_WRITE_INTERVAL = parse(Int, get(ENV, "EFM_WRITE_INTERVAL", "100"))
# EFM_SIGMA_SCALE multiplies surface tension; use 0 to disable capillarity.
EFM_SIGMA_SCALE = parse(Float64, get(ENV, "EFM_SIGMA_SCALE", "1"))
# EFM_ADAPTIVE=1 enables adaptive time stepping from the limiting Courant number.
EFM_ADAPTIVE = get(ENV, "EFM_ADAPTIVE", "1") != "0"
# EFM_MAX_CO is the adaptive target for max ordinary/film Courant number.
EFM_MAX_CO = parse(Float64, get(ENV, "EFM_MAX_CO", "0.1"))
# EFM_MIN_SHRINK is the smallest adaptive dt multiplier per iteration.
EFM_MIN_SHRINK = parse(Float64, get(ENV, "EFM_MIN_SHRINK", "0.1"))
# EFM_MAX_GROW is the largest adaptive dt multiplier per iteration.
EFM_MAX_GROW = parse(Float64, get(ENV, "EFM_MAX_GROW", "1.2"))
# XCALIBRE_EFM_DEBUG=1 prints film diagnostics every iteration.
XCALIBRE_EFM_DEBUG = get(ENV, "XCALIBRE_EFM_DEBUG", "0") != "0"
# XCALIBRE_EFM_DEBUG_INTERVAL prints diagnostics every N iterations when N > 0.
XCALIBRE_EFM_DEBUG_INTERVAL = parse(Int, get(ENV, "XCALIBRE_EFM_DEBUG_INTERVAL", "0"))
# XCALIBRE_EFM_WETTING selects hard, smooth, or allwet wetting-mask probes.
#XCALIBRE_EFM_WETTING = get(ENV, "XCALIBRE_EFM_WETTING", "hard")
XCALIBRE_EFM_WETTING = get(ENV, "XCALIBRE_EFM_WETTING", "smooth")
# XCALIBRE_EFM_WETTING_WIDTH sets the h_crit-to-width*h_crit smoothing interval.
XCALIBRE_EFM_WETTING_WIDTH = parse(Float64, get(ENV, "XCALIBRE_EFM_WETTING_WIDTH", "10"))
# XCALIBRE_EFM_FLUX_CORRECTION=0 disables pressure flux correction for A/B tests.
XCALIBRE_EFM_FLUX_CORRECTION = get(ENV, "XCALIBRE_EFM_FLUX_CORRECTION", "1") != "0"

σ *= EFM_SIGMA_SCALE

mu = 0.001003;
rho_l = 998.2; # Density of water @ 43°C kg/m3
#nu = 6.245e-7; # kinematic viscosity of water @ 43°c
nu = mu/rho_l;

cp = 4184.4 #https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html taken at 20c
T_inlet = 273+24
T_ref = 288.15

g = 9.81
gravity_vector = (g*sind(ϕ), 0.0, -g*cosd(ϕ))

gamma = 1.0
Pr = 6.9
inlet_flow_per_width = inlet_flow_rate/inlet_width # m^2/s
gravity_tangential = g*sind(ϕ)
inlet_height = (3*nu*inlet_flow_per_width/gravity_tangential)^(1/3) # m, Nusselt estimate
inlet_rate = inlet_flow_per_width/inlet_height # m/s

inlet_velocity = inlet_rate.*[1.0, 0.0, 0.0];

h_crit = 1e-10;
h_floor = 1e-15

Δt = EFM_DT # 5e-5
end;

model = Physics(
    momentum=Momentum{EFM}(;σ=σ, h_crit = h_crit, h_floor=h_floor, β=β, θm = θm,
        gravity=gravity_vector),
    time = Transient(),
    fluid = Fluid{Incompressible}(; nu = nu, rho = rho_l),
    turbulence = RANS{Laminar}(),
    energy = Energy{FilmEnergy}(Tref=T_ref),
    domain = mesh_dev
);

BCs = assign(
    region=mesh_dev,
    (
        U = [
            Dirichlet(:inlet, inlet_velocity),
            Extrapolated(:outlet),
            Wall(:inlet_sides, [0,0,0]),
            Extrapolated(:top_of_plate),
            Extrapolated(:side_1),
            Extrapolated(:side_2),
        ],
        h = [
            Dirichlet(:inlet, inlet_height),
            Extrapolated(:outlet),
            Extrapolated(:inlet_sides),
            Extrapolated(:top_of_plate),
            Extrapolated(:side_1),
            Extrapolated(:side_2)
        ],
        T = [
            #FixedTemperature(:inlet, T=T_inlet, Enthalpy(cp=cp, Tref=T_ref)),
            Dirichlet(:inlet, T_inlet),
            Zerogradient(:outlet),
            Zerogradient(:inlet_sides),
            Zerogradient(:top_of_plate),
            Zerogradient(:side_1),
            Zerogradient(:side_2)
        ]
    )
);

schemes = (
    U = Schemes(
        time=Euler,
        divergence=Upwind,
        gradient=Gauss
    ),
    h = Schemes(
        time=Euler,
        divergence=LUST
    ),
    T = Schemes(
        time=Euler,
        divergence = Upwind
    )
);

solvers = (
    U = SolverSetup(
        solver      = Bicgstab(), # Options: Gmres()
        preconditioner = Jacobi(), # Options: NormDiagonal()
        convergence = 1e-11,
        relax       = 1.0,
        rtol = 0.0,
        atol = 1e-8
    ),
    h = SolverSetup(
        solver      = Bicgstab(), # Options: Cg(), Bicgstab(), Gmres()
        preconditioner = Jacobi(), # Options: NormDiagonal()
        convergence = 1e-11,
        relax       = 1.0,
        rtol = 0.0,
        atol = 1e-8
    ),
    T = SolverSetup(
        solver      = Bicgstab(),
        preconditioner = Jacobi(),
        convergence = 1e-11,
        relax = 1.0,
        rtol = 1e-4
    )
);

adaptive = if EFM_ADAPTIVE
    AdaptiveTimeStepping(
        maxCo=EFM_MAX_CO,
        minShrink=EFM_MIN_SHRINK,
        maxGrow=EFM_MAX_GROW
    )
else
    nothing
end
begin

runtime = Runtime(
    iterations=EFM_ITERS,
    time_step=Δt,
    write_interval=EFM_WRITE_INTERVAL,
    adaptive=adaptive
)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs);

GC.gc(true)
  
initialise!(model.momentum.U, [0,0,0]);
h_init = h_floor;
initialise!(model.momentum.h, h_init)
initialise!(model.energy.Te, T_ref)

residuals = run!(model, config, inner_loops=2);

end;
