using DynamicalSystems
using LaTeXStrings, CairoMakie

# Define the dynamics
function iip_tc!(f, x, y, t)
        f[1] = x[2]*x[1]-(x[1])^2
        f[2] = 0.2
        return nothing
end

# Define the initial state and final time
x0 = [11.00, -10.0]
T = 66.5
δt = 3e-2

# Build the (non-parametric) dynamical system
tc = ContinuousDynamicalSystem(iip_tc!, x0, nothing)

# Evolve the dynamical system using the specified solver
Xt, t = trajectory(tc, T; Δt=δt)#; diffeq...)

# Print the Dataset for the time evolution
printstyled("\n", Xt, "\n"; bold=true, underline=true, color=:light_green)
#=
for j in 1:size(Xt[:],1)
        println(Xt[j])
end
=#

# Plot flow in state space
CairoMakie.activate!()
#=
fig1 = Figure(resolution = (600, 400))
ax = Axis(fig1[1, 1],
    title = "State space",
    xlabel = L"y",
    xticks = -10:5,#, [""])
    ylabel = L"x"
)
scatter!(ax, Xt[:,2], Xt[:,1], markersize=3, color = :green)
save("B-TC-11_00.png", fig1)
=#

# Create animation
points = Observable(Point2f[(x0[2],x0[1])])
fig = Figure()
ax = Axis(fig[1,1],
          title = "State space flow",
          xlabel = L"y",
          ylabel = L"x",
          xticks = -11:4, 
          yticks = -4:12
)
limits!(ax,-11,4,-4,12)
scatter!(points, markersize=3, color = :red)

frames = 1:size(Xt)[1] 
framerate = 40

record(fig, "DetTrans.gif", frames;
       framerate=framerate) do frame 
       new_point = Point2f(Xt[frame,2],Xt[frame,1])
       points[] = push!(points[], new_point)
end
