using XCALibre
# using CUDA # Uncomment to run on NVIDIA GPUs
# using AMDGPU # Uncomment to run on AMD GPUs

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS/");
#grids_dir = pkgdir(XCALibre, "Test_Meshes/");
grid = "quad.unv";
#grid = "25x25_grid.unv"
#grid = "10x10_grid.unv"
mesh_file = joinpath(grids_dir, grid);

mesh = UNV2D_mesh(mesh_file, scale=0.001);

# Select backend and setup hardware
backend = CPU();
# backend = CUDABackend() # ru non NVIDIA GPUs
# backend = ROCBackend() # run on AMD GPUs

hardware = Hardware(backend=backend, workgroup=1024);
# hardware = Hardware(backend=backend, workgroup=32) # use for GPU backends

mesh_dev = mesh; # use this line to run on CPU
# mesh_dev = adapt(backend, mesh)  # Uncomment to run on GPU 

rho_l = 100;
rho_l = 1;
rho_l = 991.07; # Density of water @ 43°C kg/m3
inlet_speed = 0.04;
h_inlet = 0.005;

velocity = inlet_speed*[1, 0.0, 0.0];
nu = 6.245e-7; # Kinematic Viscosity of water @ 43°C
Re = velocity[1]*0.01/nu;

h_crit = 1e-10;

model = Physics(
    momentum=Momentum{EFM}(; h_crit = h_crit, β=6.0, θm = 75, ϕ=90),
    time = Steady(),
    fluid = Fluid{Incompressible}(; nu = nu, rho = rho_l),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
);

BCs = assign(
    region=mesh_dev,
    (
        U = [
            Dirichlet(:inlet, velocity),
            Zerogradient(:outlet),
            Wall(:bottom, [0.0, 0.0, 0.0]),
            Wall(:top, [0.0, 0.0, 0.0])
        ],
        h = [
            Dirichlet(:inlet, h_inlet),
            Zerogradient(:outlet),
            Zerogradient(:bottom),
            Zerogradient(:top)
        ]
    )
);

schemes = (
    U = Schemes(divergence = Linear),
    h = Schemes(), # no input provided (will use defaults)
);

solvers = (
    U = SolverSetup(
        solver      = Bicgstab(), # Options: Gmres()
        preconditioner = Jacobi(), # Options: NormDiagonal()
        convergence = 1e-12,
        relax       = 0.7,
        rtol = 1e-4,
        atol = 1e-10
    ),
    h = SolverSetup(
        solver      = Bicgstab(), # Options: Cg(), Bicgstab(), Gmres()
        preconditioner = Jacobi(), # Options: NormDiagonal()
        convergence = 1e-7,
        relax       = 0.7,
        rtol = 1e-4,
        atol = 1e-10
    )
);

runtime = Runtime(iterations=2000, time_step=1, write_interval=2000)
#runtime = Runtime(iterations=20, time_step=1, write_interval=1); # hide

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs);

initialise!(model.momentum.U, velocity);
initialise!(model.momentum.h, 1);

residuals = run!(model, config);