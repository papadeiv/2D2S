using DynamicalSystems
using LaTeXStrings, CairoMakie

# Define the dynamics
function iip_bfold!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 1e-1
        return nothing
end

# Define the initial state and final time
x0 = [0.0, -1.0]
T = 14.5
δt = 5e-2

# Build the (non-parametric) dynamical system
bfold = ContinuousDynamicalSystem(iip_bfold!, x0, nothing)

# Specify solver from DifferentialEquations.jl
# diffeq = (alg=Tsit5(), reltol=1e-3, dtmax=1e-2)

# Evolve the dynamical system using the specified solver
Xt, t = trajectory(bfold, T; Δt=δt)#; diffeq...)

# Print the Dataset for the time evolution
printstyled("\n", Xt, "\n"; bold=true, underline=true, color=:light_green)

# Plot flow in state space
CairoMakie.activate!()
#=
fig1 = Figure(resolution = (600, 400))
ax = Axis(fig1[1, 1],
    title = "State space",
    xlabel = L"y",
    ylabel = L"x"
)
scatter!(ax, Xt[:,2], Xt[:,1], markersize=3, color = :green)
save("B-fold.png", fig1)
=#

# Create animation
points = Observable(Point2f[(x0[2],x0[1])])
fig = Figure()
ax = Axis(fig[1,1],
          title = "State space flow",
          xlabel = L"y",
          ylabel = L"x"
          #xticks = -1.5:0.5, 
          #yticks = -1.5:2
)
limits!(ax,-1.5,0.75,-2,1)
scatter!(points, markersize=3, color = :blue, strokecolor = :red)

frames = 1:size(Xt)[1] 
framerate = 40

record(fig, "F-bifurcation.gif", frames;
       framerate=framerate) do frame 
       new_point = Point2f(Xt[frame,2],Xt[frame,1])
       points[] = push!(points[], new_point)
end
