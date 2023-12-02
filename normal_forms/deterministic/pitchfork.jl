using DynamicalSystems
using LaTeXStrings, CairoMakie

# Define the dynamics
function iip_pf!(f, x, y, t)
        f[1] = x[2]*x[1]+y[1]*(x[1])^3
        f[2] = 0.2
        return nothing
end

# Define the initial state and final time
x0 = [0.2, -4.0]
y0 = [.5]
T = 35.0
δt = 5e-2

# Build the (non-parametric) dynamical system
pf = ContinuousDynamicalSystem(iip_pf!, x0, y0)

# Specify solver from DifferentialEquations.jl
# diffeq = (alg=Tsit5(), reltol=1e-3, dtmax=1e-2)

# Evolve the dynamical system using the specified solver
Xt, t = trajectory(pf, T; Δt=δt)#; diffeq...)

# Print the Dataset for the time evolution
printstyled("\n", Xt, "\n"; bold=true, underline=true, color=:light_green)

#=
for j in 1:size(Xt[:],1)
        println(Xt[j])
end
=#

# Plot flow in state space
#=
CairoMakie.activate!()
fig = Figure(resolution = (600, 400))
ax = Axis(fig[1, 1],
    title = "State space",
    xlabel = L"y",
    ylabel = L"x"
)
scatter!(ax, Xt[:,2], Xt[:,1], markersize=3, color = :green)
save("B-PF.png", fig)
=#

# Create animation
points = Observable(Point2f[(x0[2],x0[1])])
fig = Figure()
ax = Axis(fig[1,1],
          title = "State space flow",
          xlabel = L"y",
          ylabel = L"x",
          xticks = -5:4
)
limits!(ax,-5,4,-0.2,0.25)
scatter!(points, markersize=3, color = :blue)

frames = 1:size(Xt)[1] 
framerate = 40

record(fig, "PF-bifurcation.gif", frames;
       framerate=framerate) do frame 
       new_point = Point2f(Xt[frame,2],Xt[frame,1])
       points[] = push!(points[], new_point)
end
