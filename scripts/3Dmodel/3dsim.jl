# This is 3D simulation with unstructured mesh with Jutul package.

# Load the package and setup
import Pkg;Pkg.add("DrWatson")

using DrWatson
@quickactivate "S&P_Global_Intern"

using Jutul 
using JutulDarcy
using HYPRE
using GLMakie

nx = 100
nz = 50
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)


# Set up a 2D acquifer model
cart_dims = (nx, 1, nz)
physical_dims = (1000.0, 100.0, 50.0)
mesh = UnstructuredMesh(CartesianMesh(cart_dims, physical_dims))


points = mesh.node_points
for (i, pt) in enumerate(points)
    x, y, z = pt
    x_u = 2*π*x/1000.0
    w = 0.2
    dz = 0.05*x + w*(30*sin(2.0*x_u) + 20*sin(5.0*x_u) + 10*sin(10.0*x_u) + 5*sin(25.0*x_u))
    points[i] = pt + [0, 0, dz]
end

## 0.orginal universal permeability

# setup the simulation model
# domain = reservoir_domain(mesh, permeability = 1.0Darcy, porosity = 0.3, temperature = convert_to_si(30.0, :Celsius))
# Injector = setup_well(domain, (65, 1, 1), name = :Injector)
# model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);



## 1.create the permeability with several layers
perm1 = ones(100, 50) * 1.0Darcy
perm1[:, 1:2] *= 0.0
perm1[:, 50] *= 0.02 / 1000
perm1[:, 12] .*= 40/1000
perm1[:, 22] .*= 40/1000
perm1[:, 32] .*= 40/1000
perm1[:, 42] .*= 40/1000

poro1 = ones(100, 50) * 0.27
poro1[:, 1:2] .*= 0.0
poro1[:, 50] .*= 0.02/0.27
poro1[:, 12] .*= 0.11/0.27
poro1[:, 22] .*= 0.11/0.27
poro1[:, 32] .*= 0.11/0.27
poro1[:, 42] .*= 0.11/0.27

domain = reservoir_domain(mesh, permeability = vec(perm1), porosity = vec(poro1), temperature = convert_to_si(30.0, :Celsius))
Injector = setup_well(domain, (65, 1, 43), name = :Injector)
model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);



# ## 2.create the permeability with several layers and fault
# perm2 = ones(100, 50) * 1.0Darcy
# perm2[:, 1:2] *= 0.0
# perm2[:, 50] *= 0.02 / 1000
# perm2[:, 12] .*= 40/1000
# perm2[:, 22] .*= 40/1000
# perm2[:, 32] .*= 40/1000
# perm2[:, 42] .*= 40/1000
# perm2[50, 3:49] = 20/1000 * 1.0Darcy * ones(47)

# poro2 = ones(100, 50) * 0.27
# poro2[:, 1:2] .*= 0.0
# poro2[:, 50] .*= 0.02/0.27
# poro2[:, 12] .*= 0.11/0.27
# poro2[:, 22] .*= 0.11/0.27
# poro2[:, 32] .*= 0.11/0.27
# poro2[:, 42] .*= 0.11/0.27
# poro2[50, 3:49] = 0.05 * ones(47)

# domain = reservoir_domain(mesh, permeability = vec(perm2), porosity = vec(poro2), temperature = convert_to_si(30.0, :Celsius))
# Injector = setup_well(domain, (65, 1, 43), name = :Injector)
# model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);


# ## 3.create the permeability with several layers and baffles
# perm3 = ones(100, 50) * 1.0Darcy
# perm3[:, 1:2] *= 0.0
# perm3[:, 50] *= 0.02 / 1000
# perm3[:, 12] .*= 40/1000
# perm3[:, 22] .*= 40/1000
# perm3[:, 32] .*= 40/1000
# perm3[:, 42] .*= 40/1000

# perm3[40:45, 3:8] = 40 / 1000 * 1.0Darcy * ones(6, 6)
# perm3[55:60, 3:8] = 40 / 1000 * 1.0Darcy * ones(6, 6)
# perm3[70:75, 3:8] = 40 / 1000 * 1.0Darcy * ones(6, 6)
# perm3[80:85, 3:8] = 40 / 1000 * 1.0Darcy * ones(6, 6)

# perm3[50:53, 13:16] = 40 / 1000 * 1.0Darcy * ones(4, 4)
# perm3[60:63, 13:16] = 40 / 1000 * 1.0Darcy * ones(4, 4)
# perm3[68:71, 13:16] = 40 / 1000 * 1.0Darcy * ones(4, 4)


# perm3[20:23, 23:30] = 40 / 1000 * 1.0Darcy * ones(4, 8)
# perm3[70:73, 23:30] = 40 / 1000 * 1.0Darcy * ones(4, 8)

# perm3[30:37, 33:36] = 40 / 1000 * 1.0Darcy * ones(8, 4)
# perm3[90:97, 33:36] = 40 / 1000 * 1.0Darcy * ones(8, 4)

# perm3[45:50, 43:48] = 40 / 1000 * 1.0Darcy * ones(6, 6)
# perm3[75:80, 43:48] = 40 / 1000 * 1.0Darcy * ones(6, 6)


# poro3 = ones(100, 50) * 0.27
# poro3[:, 1:2] .*= 0.0
# poro3[:, 50] .*= 0.02/0.27
# poro3[:, 12] .*= 0.11/0.27
# poro3[:, 22] .*= 0.11/0.27
# poro3[:, 32] .*= 0.11/0.27
# poro3[:, 42] .*= 0.11/0.27

# poro3[40:45, 3:8] = 0.11 * ones(6, 6)
# poro3[55:60, 3:8] = 0.11 * ones(6, 6)
# poro3[70:75, 3:8] = 0.11 * ones(6, 6)
# poro3[80:85, 3:8] = 0.11 * ones(6, 6)

# poro3[50:53, 13:16] = 0.11 * ones(4, 4)
# poro3[60:63, 13:16] = 0.11 * ones(4, 4)

# poro3[20:23, 23:30] = 0.11 * ones(4, 8)
# poro3[70:73, 23:30] = 0.11 * ones(4, 8)

# poro3[30:37, 33:36] = 0.11 * ones(8, 4)
# poro3[90:97, 33:36] = 0.11 * ones(8, 4)

# poro3[45:50, 43:48] = 0.11 * ones(6, 6)
# poro3[75:80, 43:48] = 0.11 * ones(6, 6)


# domain = reservoir_domain(mesh, permeability = vec(perm3), porosity = vec(poro3),temperature = convert_to_si(30.0, :Celsius))
# Injector = setup_well(domain, (65, 1, 43), name = :Injector)
# model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);



# make constant pressure condition
boundary = Int[]
for cell in 1:number_of_cells(mesh)
    I, J, K = cell_ijk(mesh, cell)
    if I == 1 || I == nx
        push!(boundary, cell)
    end
end
parameters[:Reservoir][:FluidVolume][boundary] *= 1000;

# plot the model
plot_reservoir(model)


# setup schedule
nstep = 25
nstep_shut = 25
dt_inject = fill(365.0day, nstep)
pv = pore_volume(model, parameters)
inj_rate = 0.0075*sum(pv)/sum(dt_inject)

rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0],
    density = 900.0,
)


controls = Dict(:Injector => I_ctrl)
forces_inject = setup_reservoir_forces(model, control = controls)

forces_shut = setup_reservoir_forces(model)
dt_shut = fill(365.0day, nstep_shut);

dt = vcat(dt_inject, dt_shut)
forces = vcat(fill(forces_inject, nstep), fill(forces_shut, nstep_shut));

# set up initial state
state0 = setup_reservoir_state(model,
    Pressure = 200bar,
    OverallMoleFractions = [1.0, 0.0],
)

# simulate the schedule
wd, states, t = simulate_reservoir(state0, model, dt,
    parameters = parameters,
    forces = forces,
    max_timestep = 30day
)

# plot the density of the brine
using GLMakie
function plot_co2!(fig, ix, x, title = "")
    ax = Axis3(fig[ix, 1],
        zreversed = true,
        azimuth = -0.51π,
        elevation = 0.05,
        aspect = (1.0, 1.0, 0.3),
        title = title)
    plt = plot_cell_data!(ax, mesh, x, colormap = :seaborn_icefire_gradient)
    Colorbar(fig[ix, 2], plt)
end
fig = Figure(size = (900, 1200))
for (i, step) in enumerate([1, 5, nstep, nstep+nstep_shut])
    plot_co2!(fig, i, states[step][:PhaseMassDensities][1, :], "Brine density report step $step/$(nstep+nstep_shut)")
end
save("s1.png", fig)
fig

# plot the result in the interative viewer
plot_reservoir(model, states)