using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Colors
using Statistics, SpecialFunctions


#       Arnold and Boxler process 
# (stochastic transcritical bifurcation)


# Define the parameters of the process
σ = sqrt(0.8)

# Define the deterministic dynamics
function iip_deterministic!(f, x, y, t)
        f[1] = y[1]*x[1] - (x[1])^2 + (0.5)*(σ^2)*x[1]
        return nothing
end
# Define the stochastic dynamics
function iip_stochastic!(f, x, y, t)
        f[1] = +σ*x[1]
        return nothing
end

# Define initial state
x0 = [0.1]

# Define ensemble quantities
Ne = 3000

# Define parameter values
Ny = 100
counter = 1
y_values = LinRange(-0.1,+1.0,Ny)

# Define temporal evolution quantities
T = 5.00
δt = 1e-3
Nt = ceil(Int64, T/δt) + 1

# Define uninitialised array to store a single sample path 
xt = Vector{Float64}(undef, Nt)
# Define uninitialised matrix to store the ensemble entire time-series
Xt = Array{Float64, 2}(undef, Ne, Nt)
# Define uninitialised matrix to store the ensemble sample paths' final states
Et = Array{Float64, 2}(undef, Ny, Ne)
# Define uninitialised matrix to store the parametric, ensemble-averaged, sample paths
At = Array{Float64, 2}(undef, Ny, Nt)

# Loop over the parameter's values
for y in y_values
        printstyled(counter, ") Ensemble simulation for y=", y, "\n"; bold=true, underline=true, color=:light_green)
        # Perform an ensemble simulation of Ne sample paths for the current parameter value y
        for j in 1:Ne
                # Build the (non-parametric) dynamical system
                local AandB = SDEProblem(iip_deterministic!, iip_stochastic!, x0, (0.0, T), y)
                # Evolve a sample path of the stochastic process for the current parameter value 
                sol = solve(AandB, EM(), dt=δt)
                # Store the sample path into the appropriate trajectory matrix 
                global Xt[j,:] = sol[1,:]
        end
        # Compute the ensemble averaged trajectory for the current parameter value y
        for k in 1:Nt
                global xt[k] = mean(Xt[:,k]) 
        end
        # Store the ensemble's sample paths' final states into the appropriate trajectory matrix
        global Et[counter,:] = Xt[:,end]
        # Store the sample path into the appropriate trajectory matrix
        global At[counter,:] = xt 
        printstyled("   Ensemble avg=", mean(Et[counter,:]), "\n"; color=:light_blue)
        printstyled("   Temporal avg=", mean(At[counter,:]), "\n"; color=:light_red)
        # Update the counter
        global counter += 1
end

# Plot the parametrised, ensemble-averaged, sample paths 
CairoMakie.activate!()
fig1 = Figure(resolution = (1200, 900))
ax = Axis(fig1[1,1:2],
    xlabel = L"\tau",
    ylabel = L"x",
    limits = ((0,T), nothing)
)
for n in 1:Ny
        lines!(ax, LinRange(0,T,Nt), At[n,:], linewidth = 1)
end

# Plot the parametrised ensemble and temporal averaged sample paths 
ax2 = Axis(fig1[2,1:2],
    xlabel = L"y",
    ylabel = L"x",
    limits = ((y_values[1],last(y_values)), nothing)
)
# The ensemble average is computed over the final states of the sample paths of the ensemble
nsmb_avg = zeros(Ny)
# The temporal average is computed over the ensemble's averaged states across the time series
temp_avg = zeros(Ny)
for n in 1:Ny
        nsmb_avg[n] = mean(Et[n,:])
        temp_avg[n] = mean(At[n,:])
end
lines!(ax2, y_values, nsmb_avg, color = :red, linewidth = 1)
lines!(ax2, y_values, temp_avg, color = :blue, linewidth = 1)

# Plot the (deterministic) critical manifold
y_pos = range(1,0,length=10)
y_neg = range(-1,0,length=10)
stable_pos = 0.0 .+ 1.0*y_pos 
stable_neg = 0.0 .+ 0.0*y_neg
unstable_pos = 0.0 .+ 0.0*y_pos 
unstable_neg = 0.0 .+ 1.0*y_neg
lines!(y_pos, stable_pos, color = :black, linewidth = 0.5)
lines!(y_pos, unstable_pos, color = :black, linewidth = 0.5, linestyle = :dash)
lines!(y_neg, stable_neg, color = :black, linewidth = 0.5)
lines!(y_neg, unstable_neg, color = :black, linewidth = 0.5, linestyle = :dash)

# Plot P-bifurcation point of the transcritical normal form
yP = (σ^2)/2.0
points = Observable(Point2f[(-yP, 0.0)])
new_point = Point2f(yP, 0.0)
points[] = push!(points[], new_point)
scatter!(ax2, points, markersize=10, color = :blue)
text!(L"y^{+}_{P}", position = (yP,-0.25), align = (:center,:center), color = :blue, fontsize = 20)
text!(L"y^{-}_{P}", position = (-yP,-0.25), align = (:center,:center), color = :blue, fontsize = 20)
sing = Observable(Point2f[(0.0,0.0)])
scatter!(ax2, sing, markersize=10, color = :black)
text!(L"y_{D}", position = (0.0,-0.25), align = (:center,:center), color = :black, fontsize = 20)

# Plot distributions of parametrised ensemble sample paths at the last state 
#=
ax3 = Axis(fig1[3,1],
         title = L"y=-0.8",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau=5.0}",
         yticks = [0,0.035,0.07],
         limits = (0,1,0,0.075)
)
hist!(ax3, Et[9,:], bins = 100, normalization = :probability, color = :red, strokecolor = :black, strokewidth = 1)

ax4 = Axis(fig1[3,2],
         title = L"y=-0.2",
         xlabel = L"x_{\tau}",
         yticks = [0,0.035,0.07],
         limits = (0,1,0,0.075)
)
hist!(ax4, Et[32,:], bins = 100, normalization = :probability, color = :red, strokecolor = :black, strokewidth = 1)
=#

ax5 = Axis(fig1[3,1],
         title = L"y=+0.2",
         limits = (0,1,nothing,nothing)
)
hist!(ax5, Et[28,:], bins = 50, color = :red, strokecolor = :black, strokewidth = 1)

ax6 = Axis(fig1[3,2],
         title = L"y=+0.8",
         limits = (0,2,nothing,nothing)
)
hist!(ax6, Et[82,:], bins = 50, color = :red, strokecolor = :black, strokewidth = 1)

# Plot the theoretical and numerical variance of the parametrised distribution 
fig2 = Figure(resolution = (1200, 900))
ax1 = Axis(fig2[1,1],
           xlabel = L"y",
           ylabel = L"v_s",
           limits = ((-1.0,+1.0), (0,2))
)
# Numerical values computed from time-series
Var = zeros(Float64, Ny)
for n in 1:Ny
        Var[n] = var(Et[n,:])
end
# Theoretical values [Kuehn, 2011] 
#Var_exact = time.^2 .+ 0.5*(σ^2).*time .- 4 .^(.-(1/σ^2).*time).*(1/σ^2).^(-(2/σ^2).*time).*(σ^2).*gamma.(1 .+(2/σ^2).*time) .+ 2 .^(.-(4/σ^2).*time.^2).*(1/σ^2).^(-(4/σ^2).*time).*(gamma.((2/σ^2).*time)).^2
l1 = lines!(ax1, y_values, Var, color = :purple)
#l2 = lines!(ax3, time, Var_exact, color = :blue, linestyle = :dash)

# Export the figures
save("SteadyStateFokkerPlanck.png", fig1)
save("Variance.png", fig2)
