using DynamicalSystems
using LaTeXStrings, CairoMakie

# LAYER EQUATION

# Define the dynamics
function iip_bfold!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.0 
        return nothing
end

# Define the initial state and final time
x0 = [-0.50, -1.00]
x1 = [1.2, -0.875]
x2 = [0.00, -0.75]
x3 = [1.125, -0.625]
x4 = [0.25, -0.50]
x5 = [1.00, -0.375]
T = 10.0
δt = 1e-1

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_bfold!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X0 = Xt[:,1]
Y0 = Xt[:,2]
bfold = ContinuousDynamicalSystem(iip_bfold!, x1, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X1 = Xt[:,1]
Y1 = Xt[:,2]
bfold = ContinuousDynamicalSystem(iip_bfold!, x2, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X2 = Xt[:,1]
Y2 = Xt[:,2]
bfold = ContinuousDynamicalSystem(iip_bfold!, x3, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X3 = Xt[:,1]
Y3 = Xt[:,2]
bfold = ContinuousDynamicalSystem(iip_bfold!, x4, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X4 = Xt[:,1]
Y4 = Xt[:,2]
bfold = ContinuousDynamicalSystem(iip_bfold!, x5, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X5 = Xt[:,1]
Y5 = Xt[:,2]
# Simulation parameters
Nt = size(X0,1)
y_values = LinRange(-1.125,0.00,200)

# Plot fast flow in state space
CairoMakie.activate!()
fig1 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax = Axis(fig1[1, 1],
    backgroundcolor = :transparent,
    xlabel = L"y",
    ylabel = L"x",
    limits = ((y_values[1],y_values[end]+0.05), (-1.25,+1.25))
)
lines!(ax, Y0, X0, color = :blue, linewidth = 1.0)
lines!(ax, Y1, X1, color = :blue, linewidth = 1.0)
lines!(ax, Y2, X2, color = :blue, linewidth = 1.0)
lines!(ax, Y3, X3, color = :blue, linewidth = 1.0)
lines!(ax, Y4, X4, color = :blue, linewidth = 1.0)
lines!(ax, Y5, X5, color = :blue, linewidth = 1.0)
# Plot the (deterministic) critical manifold
stable = sqrt.(-y_values)
unstable = -sqrt.(-y_values)
lines!(ax, y_values, stable, color = :black, linewidth = 1.5)
lines!(ax, y_values, unstable, color = :black, linewidth = 1.5, linestyle = :dash)
# Plot ICs of each trajectory
scatter!(ax, x0[2], x0[1], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
#text!(x0[2], x0[1], text="(".*string.(round(x0[2], digits=2)).*",".* string.(round(x0[1], digits=2)).*")", align = (:center, :center))
scatter!(ax, x1[2], x1[1], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, x2[2], x2[1], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, x3[2], x3[1], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, x4[2], x4[1], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, x5[2], x5[1], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot the final states of each trajectory
scatter!(ax, Y0[end], X0[end], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, Y1[end], X1[end], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, Y2[end], X2[end], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, Y3[end], X3[end], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, Y4[end], X4[end], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, Y5[end], X5[end], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)

# Export the results
#save("../results/F-bifurcation.png", fig1)

# FAST-SLOW SYSTEM

# Define the dynamics
function iip_fs!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.01 
        return nothing
end

# Define the initial state and final time
x0 = [1.00,-1.00]
T = 70.0
δt = 1e0

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_fs!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X = Xt[:,1]
Y = Xt[:,2]

# Plot fast flow in state space
CairoMakie.activate!()
fig2 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax1 = Axis(fig2[1, 1],
    backgroundcolor = :transparent,
    xlabel = L"y",
    ylabel = L"x",
    limits = ((y_values[1],y_values[end]+0.05), (-1.25,+1.25))
)
# Plot the (deterministic) critical manifold
stable = sqrt.(-y_values)
unstable = -sqrt.(-y_values)
lines!(ax1, y_values, stable, color = :black, linewidth = 1.5)
lines!(ax1, y_values, unstable, color = :black, linewidth = 1.5, linestyle = :dash)
# Plot the flow in state space
scatter!(ax1, Y, X, color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax1, Y0, X0, color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot the end states of each trajectory
scatter!(ax1, Y0[1], X0[1], color = :purple, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax1, Y[1], X[1], color = :yellow, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax1, Y[end], X[end], color = :green, strokecolor = :black, strokewidth = 1.5, markersize = 9)

# Export the results
save("../results/F-bifurcation.png", fig2)

# Create animation
#=
points = Observable(Point2f[(x0[2],x0[1])])
fig = Figure(; size = (600,400))
ax = Axis(fig[1,1],
          title = "State space flow",
          xlabel = L"y",
          ylabel = L"x"
)
limits!(ax,-1.5,0.75,-2,1)
scatter!(points, markersize=3, color = :blue, strokecolor = :red)

frames = 1:size(Xt)[1] 
framerate = 40

record(fig, "../results/DetFold.gif", frames;
       framerate=framerate) do frame 
       new_point = Point2f(Xt[frame,2],Xt[frame,1])
       points[] = push!(points[], new_point)
end
=#
