using DynamicalSystems
using LaTeXStrings, CairoMakie

# Define the dynamics
function iip_bfold!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.02
        return nothing
end

# Define the initial state and final time
x0 = [1.2, -0.6]
T = 37.75
δt = 5e-2

# Build the (non-parametric) dynamical system
bfold = ContinuousDynamicalSystem(iip_bfold!, x0, nothing)

# Specify solver from DifferentialEquations.jl
# diffeq = (alg=Tsit5(), reltol=1e-3, dtmax=1e-2)

# Evolve the dynamical system using the specified solver
Xt, t = trajectory(bfold, T; Δt=δt)#; diffeq...)

# Print the Dataset for the time evolution
printstyled("\n", Xt, "\n"; bold=true, underline=true, color=:light_green)
#=
for j in 1:size(Xt[:],1) println(Xt[j])
end
=#

# Plot flow in state space
CairoMakie.activate!()
fig1 = Figure(resolution = (1000, 800))
ax = Axis(fig1[1, 1],
    title = "State space",
    xlabel = L"y",
    ylabel = L"x"
)
limits!(ax,-0.7,0.2,-1.2,1.2)
scatter!(ax, Xt[:,2], Xt[:,1], markersize=3, color = :green)
y_critical = range(-1,0,length=1000)
x_attr_critical = sqrt.(-y_critical)
lines!(ax, y_critical, x_attr_critical, color = :black, linewidth = 1, linestyle = :solid)
x_repl_critical = -sqrt.(-y_critical)
lines!(ax, y_critical, x_repl_critical, color = :black, linewidth = 1, linestyle = :dash)
save("B-fold.png", fig1)

# Create animation
#=
points = Observable(Point2f[(x0[2],x0[1])])
fig = Figure()
ax = Axis(fig[1,1],
          title = "State space flow",
          xlabel = L"y",
          ylabel = L"x"
          #xticks = -11:1, 
          #yticks = -3:3
)
limits!(ax,-0.7,0.2,-1.2,1.2)
scatter!(points, markersize=4, color = :red, strokecolor = :black)
y_critical = range(-1,0,length=1000)
x_attr_critical = sqrt.(-y_critical)
lines!(ax, y_critical, x_attr_critical, color = :black, linewidth = 1, linestyle = :solid)
x_repl_critical = -sqrt.(-y_critical)
lines!(ax, y_critical, x_repl_critical, color = :black, linewidth = 1, linestyle = :dash)

frames = 1:size(Xt)[1] 
framerate = 40

record(fig, "Example 2.3.gif", frames;
       framerate=framerate) do frame 
       new_point = Point2f(Xt[frame,2],Xt[frame,1])
       points[] = push!(points[], new_point)
end
=#
