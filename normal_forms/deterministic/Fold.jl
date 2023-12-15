using DynamicalSystems, Colors
using LaTeXStrings, CairoMakie

# FAST-SLOW EQUATION

# Define the dynamics
function iip_bfold!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.01
        return nothing
end

# Define the initial state and final time
x0 = [-0.50, -1.00]
T = 110.10
δt = 1e-1
y_values = LinRange(-1.125,0.00,200)

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_bfold!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X0 = Xt[:,1]
Y0 = Xt[:,2]

# Plot the (deterministic) critical manifold
CairoMakie.activate!()
fig1 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax = Axis(fig1[1, 1],
    backgroundcolor = :transparent,
    xlabel = L"y",
    ylabel = L"x",
    limits = ((y_values[1],y_values[end]+0.15), (-1.25,+1.25))
)
stable = sqrt.(-y_values)
unstable = -sqrt.(-y_values)
lines!(ax, y_values, stable, color = :black, linewidth = 1.5)
lines!(ax, y_values, unstable, color = :black, linewidth = 1.5, linestyle = :dash)

# SINGULAR LIMIT - Trajectory 1

# Define the dynamics
function iip_fs1!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.0 
        return nothing
end

# Define the initial state and final time
x0 = [-0.50,-1.00]
T = 10.0
δt = 1e-1

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_fs1!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X1 = Xt[:,1]
Y1 = Xt[:,2]

# SINGULAR LIMIT - Trajectory 2

# Define the dynamics
function iip_fs2!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.0001 
        return nothing
end

# Define the initial state and final time
x0 = [+1.00,-1.00]
T = 10000.0
δt = 100e0

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_fs2!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X2 = Xt[:,1]
Y2 = Xt[:,2]

# SINGULAR LIMIT - Trajectory 3

# Define the dynamics
function iip_fs3!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.0 
        return nothing
end

# Define the initial state and final time
x0 = [-0.01,0.0]
T = 100.0
δt = 1e-1

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_fs3!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X3 = Xt[:,1]
Y3 = Xt[:,2]

# Plot singular limit trajectory 1 
scatter!(ax, Y1, X1, color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot singular limit trajectory 2 
scatter!(ax, Y2, X2, color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot singular limit trajectory 3 
scatter!(ax, Y3, X3, color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot the fast-slow trajectory 
lines!(ax, Y0, X0, color = :orange, linewidth = 2.5)
# Export the results
save("../results/F-bifurcation.png", fig1)
