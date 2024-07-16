# This is a toy case for exploring the CO2 plume behavior 
# on different geological structure 

using JutulDarcyRules
using CSV 
using PyPlot
using PyCall
using DataFrames

# enable making videos
@pyimport matplotlib.animation as anim 

# Set up the parameters for simulation
ϕ = Matrix(CSV.read("poro.csv", DataFrame; header=false))
ϕ = transpose(ϕ)
K = Matrix(CSV.read("perm.csv", DataFrame; header=false))
K = transpose(K)
n = (18, 1, 18)
d = (20.0, 1000.0, 20.0)
h = 1000.0

# hydraulic pressure
P_hydrau = (repeat(collect(1:18), 1, 18) * d[3] .+ h) * JutulDarcyRules.ρH2O * 9.807   

# injection location
inj_loc_grid = (10, 1, 13)
inj_loc = (inj_loc_grid[1], 1, inj_loc_grid[3]) .* d

# injection rate
irate = [0.05]

# plot for permeability and porosity
PyPlot.rc("font", family="serif")
PyPlot.rc("xtick", labelsize=15)
PyPlot.rc("ytick", labelsize=15)
fig=figure(figsize=(8, 8))
imshow(K', extent=(0, (n[1]-1)*d[1], h+(n[3]-1)*d[3], h), aspect="equal", cmap="Blues")
im_ratio = size(K, 2)*d[3] / (size(K, 1)*d[1])
clb = colorbar(fraction=0.046*im_ratio, pad=0.04)
clb[:ax][:set_title]("Md", fontsize=15)

# label x and y axis
xlabel("X[m]", fontsize=15)
ylabel("Depth[m]", fontsize=15)
title("Permeability", fontsize=20)

# plot injection location
plot(inj_loc[1]-d[1], inj_loc[3]+h, "xk", markersize=15)
plot(inj_loc[1]-d[1], inj_loc[3]+h+d[3], "xk", markersize=15)
plot(inj_loc[1]-d[1], inj_loc[3]+h+2*d[3], "xk", markersize=15)
text(inj_loc[1]+d[1]-d[1], inj_loc[3]+h+d[3], c="k", "Injection", fontsize=20)
plt.tight_layout(pad=0.6, w_pad=0.5, h_pad=1.0)

fig.savefig("perm")
clf()


fig=figure(figsize=(8, 8))
imshow(ϕ', extent=(0, (n[1]-1)*d[1], h+(n[3]-1)*d[3], h), cmap="Blues")
im_ratio = size(ϕ, 2)*d[3] / (size(ϕ, 1)*d[1])
clb = colorbar(fraction=0.046*im_ratio, pad=0.04)

# label x and y axis
xlabel("X[m]", fontsize=15)
ylabel("Depth[m]", fontsize=15)
title("Porosity", fontsize=20)

# plot injection location
plot(inj_loc[1]-d[1], inj_loc[3]+h, "xk", markersize=15)
plot(inj_loc[1]-d[1], inj_loc[3]+h+d[3], "xk", markersize=15)
plot(inj_loc[1]-d[1], inj_loc[3]+h+2*d[3], "xk", markersize=15)
text(inj_loc[1]+d[1]-d[1], inj_loc[3]+h+d[3], c="k", "Injection", fontsize=20)
plt.tight_layout(pad=0.6, w_pad=0.5, h_pad=1.0)

fig.savefig("poro")
clf()

# plot initial hydraulic pressure
fig = figure(figsize=(8, 8))
imshow(P_hydrau / 10^6, extent=(0, (n[1]-1)*d[1], h+(n[3]-1)*d[3], h), cmap="Blues")
im_ratio = size(K, 2)*d[3] / (size(K, 1)*d[1])
clb = colorbar(fraction=0.046*im_ratio, pad=0.04)
clb[:ax][:set_title]("MPa", fontsize=15)

# label x and y axis
xlabel("X[m]", fontsize=15)
ylabel("Depth[m]", fontsize=15)
title("Hydraulic Pressure", fontsize=20)

# plot injection location
plot(inj_loc[1]-d[1], inj_loc[3]+h, "xk", markersize=15)
plot(inj_loc[1]-d[1], inj_loc[3]+h+d[3], "xk", markersize=15)
plot(inj_loc[1]-d[1], inj_loc[3]+h+2*d[3], "xk", markersize=15)
text(inj_loc[1]+d[1]-d[1], inj_loc[3]+h+d[3], c="k", "Injection", fontsize=20)
plt.tight_layout(pad=0.6, w_pad=0.5, h_pad=1.0)

fig.savefig("hydraupress")
clf()


# start the reservior simulation for the above parameter setting

# convert permeability to the millidarcy unit
K = K * JutulDarcyRules.md

# simulate for 10 years
tstep = 365 * ones(10)

# model = jutulModel(n, d, vec(Float64.(ϕ)), K1to3(K; kvoverkh=0.36), h)
model = jutulModel(n, d, 0.25, K1to3(K; kvoverkh=0.36); h=h)
f = jutulVWell(irate[1], [(inj_loc[1], inj_loc[2])]; startz=[inj_loc[3]], endz=[inj_loc[3]+3*d[3]])
# f = jutulSource(irate[1], inj_loc)
S = jutulModeling(model, tstep)
Trans = KtoTrans(CartesianMesh(model), K1to3(K; kvoverkh=0.36))
@time states = S(log.(Trans), f)                      

# get the time-varying saturation and pressure
Sat = [reshape(states.states[i][1:n[1]*n[3]], n[1], n[end]) for i in 1:10]
P = [reshape(states.states[i][n[1]*n[3]+1:end], n[1], n[end]) for i in 1:10]


# plot the saturation and pressure
for i in 1:10
    figure(figsize=(9, 16))
    subplot(2,1,1)
    sat = Sat[i]
    imshow(sat', extent=(0, (n[1]-1)*d[1], h+(n[3]-1)*d[3], h), vmin=0, vmax=1, cmap="Blues")
    clb = colorbar(fraction=0.046*im_ratio, pad=0.04)

    plot(inj_loc[1]-d[1], inj_loc[3]+h, "xk", markersize=15)
    plot(inj_loc[1]-d[1], inj_loc[3]+h+d[3], "xk", markersize=15)
    plot(inj_loc[1]-d[1], inj_loc[3]+h+2*d[3], "xk", markersize=15)
    text(inj_loc[1]+d[1]-d[1], inj_loc[3]+h+d[3], c="k", "Injection", fontsize=20)

    # label x and y axis
    xlabel("X[m]", fontsize=15)
    ylabel("Depth[m]", fontsize=15)
    title("Saturation at year " * string(i), fontsize=20, y=1.05)

    subplot(2,1,2)
    p = P[i] - P_hydrau'
    imshow(p'/10^6, extent=(0, (n[1]-1)*d[1], h+(n[3]-1)*d[3], h), cmap="Blues")
    clb = colorbar(fraction=0.046*im_ratio, pad=0.04)
    clb[:ax][:set_title]("MPa", fontsize=15)

    plot(inj_loc[1]-d[1], inj_loc[3]+h, "xk", markersize=15)
    plot(inj_loc[1]-d[1], inj_loc[3]+h+d[3], "xk", markersize=15)
    plot(inj_loc[1]-d[1], inj_loc[3]+h+2*d[3], "xk", markersize=15)
    text(inj_loc[1]+d[1]-d[1], inj_loc[3]+h+d[3], c="k", "Injection", fontsize=20)

    # label x and y axis
    xlabel("X[m]", fontsize=15)
    ylabel("Depth[m]", fontsize=15)
    title("Pressure difference between current and hydraulic at year " * string(i), fontsize=20, y=1.05)

    plt.tight_layout(pad=0.6, w_pad=0.5, h_pad=1.0)
    savefig("sat_press_year" * string(i))

    clf()

end



# make videos for simulation result of saturation and pressure 

# 10 figures
frames = 9

fig, axes = subplots(2, 1, figsize=(9, 16))

ax1, ax2 = axes

sat = Sat[1]
im_sat = imshow(sat', extent=(0, (n[1]-1)*d[1], h+(n[3]-1)*d[3], h), vmin=0, vmax=1, cmap="Blues")

ax1.plot(inj_loc[1]-d[1], inj_loc[3]+h-d[3], "xk", markersize=15)
ax1.plot(inj_loc[1]-d[1], inj_loc[3]+h+d[3]-d[3], "xk", markersize=15)
ax1.plot(inj_loc[1]-d[1], inj_loc[3]+h+2*d[3]-d[3], "xk", markersize=15)
ax1.text(inj_loc[1]+d[1]-d[1], inj_loc[3]+h+d[3]-d[3], c="k", "Injection", fontsize=20)

fig.colorbar(im_sat, fraction=0.046*im_ratio, pad=0.04)
ax1.set_xlabel("X[m]", fontsize=15)
ax1.set_ylabel("Depth[m]", fontsize=15)
ax1.set_title("Saturation at year " * string(1), fontsize=20, y=1.05)


p = P[1] - P_hydrau'
im_p = ax2.imshow(p'/10^6, extent=(0, (n[1]-1)*d[1], h+(n[3]-1)*d[3], h), cmap="Blues")
clb = colorbar(fraction=0.046*im_ratio, pad=0.04)
clb[:ax][:set_title]("MPa", fontsize=15)

ax2.plot(inj_loc[1]-d[1], inj_loc[3]+h-d[3], "xk", markersize=15)
ax2.plot(inj_loc[1]-d[1], inj_loc[3]+h+d[3]-d[3], "xk", markersize=15)
ax2.plot(inj_loc[1]-d[1], inj_loc[3]+h+2*d[3]-d[3], "xk", markersize=15)
ax2.text(inj_loc[1]+d[1]-d[1], inj_loc[3]+h+d[3]-d[3], c="k", "Injection", fontsize=20)

# label x and y axis
ax2.set_xlabel("X[m]", fontsize=15)
ax2.set_ylabel("Depth[m]", fontsize=15)
ax2.set_title("Pressure difference between current and hydraulic at year " * string(1), fontsize=20, y=1.05)
plt.tight_layout(pad=0.6, w_pad=0.5, h_pad=1.0)


# function to update the frame
function update(frame_number)
    sat = Sat[frame_number+1]
    im_sat.set_array(sat)

    p = P[frame_number+1] - P_hydrau'
    im_p.set_array(p)

    return (im_sat, im_p)
end


# create animation object
ani = anim.FuncAnimation(fig, update, frames=9, interval=200, blit=true)

# save the animation
ani.save("sat_press", writer="ffmpeg")

# close the figure
close(fig)

