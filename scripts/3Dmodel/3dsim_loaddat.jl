# This is 3D simulation with unstructured mesh with Jutul package.

# Load the package and setup
import Pkg;Pkg.add("DrWatson")

using DrWatson
@quickactivate "S&P_Global_Intern"

using Jutul 
using JutulDarcy
using HYPRE
using GLMakie
using CSV
using DataFrames


nx = 31
ny = 30
nz = 10
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)

# Set up a 2D acquifer model
cart_dims = (nx, ny, nz)
physical_dims = (1000.0, 1000.0, 1000.0)
mesh = UnstructuredMesh(CartesianMesh(cart_dims, physical_dims))

# Read the xyz position from 
top = Matrix(CSV.read("../../data/CCS2_small_xyz.dat", DataFrame; header=false))
base = Matrix(CSV.read("../../data/CCS2_small_base_xyz.dat", DataFrame; header=false))

# using Geodesy

# Gulf of mexico is in UTM zone 15, north hemispehre
# utm_zone = UTMZ(15, true, wgs84)



# depth = base[1, 3] - base[(nx + 1) * (ny + 1), 3]
# # using Ranges
# depth_arr = LinRange(0.0, depth, 11) .- base[1, 3]


# depth_pair = zeros(2, (nx + 1) * (ny + 1))
# depth_pair[1, 1:(nx + 1) * (ny + 1)] = (-1) * top[1:(nx + 1) * (ny + 1), 3]
# depth_pair[2, 1:(nx + 1) * (ny + 1)] = (-1) * base[1:(nx + 1) * (ny + 1), 3]



points = mesh.node_points
for (i, pt) in enumerate(points)
    # x, y, z = pt
    # x_u = 2*π*x/1000.0
    # w = 0.2
    # dz = 0.05*x + w*(30*sin(2.0*x_u) + 20*sin(5.0*x_u) + 10*sin(10.0*x_u) + 5*sin(25.0*x_u))
    # points[i] = pt + [0, 0, dz]

    # xy_utm = base[i % ((nx + 1) * (ny + 1)) == 0 ? (nx + 1) * (ny + 1) : i % ((nx + 1) * (ny + 1)), :]
    # x_meters, y_meters = utm_to_meters(utm_zone, xy_utm[1], xy_utm[2])
    # points[i] = [x_meters, y_meters, depth_arr[i == (nx + 1) * (ny + 1) * (nz + 1)  ? nz + 1 : div(i, (nx + 1) * (ny + 1)) + 1]]
    
    xyz_top = top[i % ((nx + 1) * (ny + 1)) == 0 ? (nx + 1) * (ny + 1) : i % ((nx + 1) * (ny + 1)), :]
    xyz_base = base[i % ((nx + 1) * (ny + 1)) == 0 ? (nx + 1) * (ny + 1) : i % ((nx + 1) * (ny + 1)), :]
    rank = div(i, (nx + 1) * (ny + 1)) == 11 ? 11 : div(i, (nx + 1) * (ny + 1)) + 1 
    depth_arr = LinRange((-1) * xyz_top[3], (-1) * xyz_base[3], 11) 
    # points[i] = [xyz_base[1] - base[1, 1], xyz_base[2] -  base[(nx + 1) * (ny + 1), 2], depth_arr[i == (nx + 1) * (ny + 1) * (nz + 1)  ? nz + 1 : div(i, (nx + 1) * (ny + 1)) + 1]]
    points[i] = [xyz_base[1] - base[1, 1], xyz_base[2] -  base[(nx + 1) * (ny + 1), 2], depth_arr[rank]]

end



# 0.orginal universal permeability

# setup the simulation model
domain = reservoir_domain(mesh, permeability = 1.0Darcy, porosity = 0.3, temperature = convert_to_si(30.0, :Celsius))

Injector_xy = [302600 - base[1, 1], 3176650 - base[(nx + 1) * (ny + 1), 2]]

Injector_cor = (round((Injector_xy[1] - points[1][1]) / (points[nx + 2][1] - points[1][1])), 
round((points[1][2] - Injector_xy[2]) / (points[1][2] - points[2][2])))


# Injector = setup_well(domain, (65, 1, 1), name = :Injector)
Injector = setup_well(domain, (Int(Injector_cor[1]), Int(Injector_cor[2]), nz - div(nz, 4)), name = :Injector)
model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);





# ## 1.create the permeability with several layers
# perm1 = ones(nx, ny, nz) * 1.0Darcy
# perm1[:, :, 1] *= 0.0
# perm1[:, :, nz] *= 0.02 / 1000
# perm1[:, :, 3] .*= 40/1000
# perm1[:, :, 5] .*= 40/1000
# perm1[:, :, 7] .*= 40/1000


# poro1 = ones(nx, ny, nz) * 0.27
# poro1[:, :, 1] .*= 0.0
# poro1[:, :, nz] .*= 0.02/0.27
# poro1[:, :, 3] .*= 0.11/0.27
# poro1[:, :, 5] .*= 0.11/0.27
# poro1[:, :, 7] .*= 0.11/0.27

# domain = reservoir_domain(mesh, permeability = vec(perm1), porosity = vec(poro1), temperature = convert_to_si(30.0, :Celsius))

# Injector_xy = [302600 - base[1, 1], 3176650 - base[(nx + 1) * (ny + 1), 2]]

# Injector_cor = (round((Injector_xy[1] - points[1][1]) / (points[nx + 2][1] - points[1][1])), 
# round((points[1][2] - Injector_xy[2]) / (points[1][2] - points[2][2])))


# # Injector = setup_well(domain, (65, 1, 1), name = :Injector)
# Injector = setup_well(domain, (Int(Injector_cor[1]), Int(Injector_cor[2]), nz - div(nz, 4)), name = :Injector)
# model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);






# ## 2.create the permeability with several layers and fault
# perm2 = ones(nx, ny, nz) * 1.0Darcy
# perm2[:, :, 1] *= 0.0
# perm2[:, :, nz] *= 0.02 / 1000
# perm2[:, :, 3] .*= 40/1000
# perm2[:, :, 5] .*= 40/1000
# perm2[:, :, 7] .*= 40/1000


# Injector_xy = [302600 - base[1, 1], 3176650 - base[(nx + 1) * (ny + 1), 2]]

# Injector_cor = (round((Injector_xy[1] - points[1][1]) / (points[nx + 2][1] - points[1][1])), 
# round((points[1][2] - Injector_xy[2]) / (points[1][2] - points[2][2])))


# perm2[Int(Injector_cor[1]) - 1:Int(Injector_cor[1]) + 2, Int(Injector_cor[2]) - 1, 2:10] = 20/1000 * 1.0Darcy * ones(4, 9)
# perm2[Int(Injector_cor[1]) - 1, Int(Injector_cor[2]) - 1:Int(Injector_cor[2]) + 2, 2:10] = 20/1000 * 1.0Darcy * ones(4, 9)

# poro2 = ones(nx, ny, nz) * 0.27
# poro2[:, :, 1] .*= 0.0
# poro2[:, :, nz] .*= 0.02/0.27
# poro2[:, :, 3] .*= 0.11/0.27
# poro2[:, :, 5] .*= 0.11/0.27
# poro2[:, :, 7] .*= 0.11/0.27

# poro2[Int(Injector_cor[1]) - 1:Int(Injector_cor[1]) + 2,  Int(Injector_cor[2]) - 1 , 2:10] = 0.05 * ones(4, 9)
# poro2[Int(Injector_cor[1]) - 1,  Int(Injector_cor[2]) - 1:Int(Injector_cor[2]) + 2, 2:10] = 0.05 * ones(4, 9)

# domain = reservoir_domain(mesh, permeability = vec(perm2), porosity = vec(poro2), temperature = convert_to_si(30.0, :Celsius))
# Injector = setup_well(domain, (Int(Injector_cor[1]), Int(Injector_cor[2]), nz - div(nz, 4)), name = :Injector)
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




# # # modify this, this way to setting up pressure condition does not work in real 3D

# # make constant pressure condition
# boundary = Int[]
# for cell in 1:number_of_cells(mesh)
#     I, J, K = cell_ijk(mesh, cell)
#     if I == 1 || I == nx
#         push!(boundary, cell)
#     end
# end
# parameters[:Reservoir][:FluidVolume][boundary] *= 1000;


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
save("loads0.png", fig)
fig

# plot the result in the interative viewer
plot_reservoir(model, states)