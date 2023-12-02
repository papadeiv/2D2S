using DynamicalSystems
using LaTeXStrings, CairoMakie

σ = sqrt(0.8)

# Define the dynamics
function iip_FP!(f, p, y, t)
        f[1] = p[2] 
        f[2] = ()*p[1] + 2*(y[1]/(σ^2)*t^2) - 4/(σ^2)*t + 1/t^2)*p[2]
        return nothing
end

# Define the initial state and final time
x0 = [0.5, 1.0]
T = 1.0
δt = 1e-3

# Build the (non-parametric) dynamical system
FP = ContinuousDynamicalSystem(iip_FP!, x0, nothing)

# Evolve the dynamical system using the specified solver
Xt, t = trajectory(FP, T; Δt=δt)#; diffeq...)

# Print the Dataset for the time evolution
printstyled("\n", Xt, "\n"; bold=true, underline=true, color=:light_green)

# Plot flow in state space
CairoMakie.activate!()
fig1 = Figure(resolution = (1200, 800))
ax = Axis(fig1[1, 1],
    title = "State space",
    xlabel = L"y",
    ylabel = L"x"
)
scatter!(ax, Xt[:,1], Xt[:,2], markersize=3, color = :green)
save("FokkerPlanck.png", fig1)
