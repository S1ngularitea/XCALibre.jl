using XCALibre
#https://essay.utwente.nl/fileshare/file/63314/BSc_report_Peter_Puttkammer.pdf
grids_dir = pkgdir(XCALibre, "test_meshes/")
grid = "Mesh_2.unv"
mesh_file = joinpath(grids_dir, grid)

mesh = UNV2D_mesh(mesh_file, scale=0.001)

backend = CPU()

hardware = Hardware(backend=backend, workgroup=1024)

mesh_dev = mesh

# define flow conditions standins used for now
rho = 1.2041 # kg/m^3
mu = 1.827e-5 # Pas
nu = mu/rho
free_stream_velocity = 10
#free_stream_velocity = 20
velocity = [free_stream_velocity,0,0]
pressure = 0
x = 0.095
Re = (rho*free_stream_velocity*x)/mu

model = Physics(
    time = Steady(),
    fluid = Fluid{Incompressible}(nu = nu),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
)

BCs = assign(
    region=mesh_dev,
    (
        U = [
            Dirichlet(:Inlet, velocity),
            Extrapolated(:Outlet),
            Wall(:Wall, [0.0, 0.0, 0.0]),
            Extrapolated(:Top)
        ],
        p = [
            Extrapolated(:Inlet), # check this
            Dirichlet(:Outlet, pressure), # and this
            Wall(:Wall),
            Extrapolated(:Top)
        ]
    )
)

schemes = (
    U = Schemes(divergence = Linear),
    p = Schemes()
)

solvers = (
    U = SolverSetup(
        solver = Bicgstab(),
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax = 0.7,
        rtol = 1e-1
    ),
    p = SolverSetup(
        solver = Cg(),
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax = 0.7,
        rtol = 1e-1
    )
)

runtime = Runtime(iterations=2000, time_step = 1, write_interval=200)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs
)

initialise!(model.momentum.U, velocity)
initialise!(model.momentum.p, 0.0)

residuals = run!(model, config)