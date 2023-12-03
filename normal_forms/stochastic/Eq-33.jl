using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Colors
using Statistics, SpecialFunctions, ProgressMeter

#                  Arnold and Boxler process 
# (stochastic transcritical bifurcation with multiplicative noise)

# Define the parameters of the process
σ = sqrt(0.8)

# Define initial state
x0 = [0.1]

# Define temporal evolution quantities
T = 10.00
δt = 1e-3

# Define the deterministic dynamics
function iip_det!(f, x, y, t)
        f[1] = y[1]*x[1] - (x[1])^2 + (0.5)*(σ^2)*x[1]
        return nothing
end
# Define the stochastic dynamics
function iip_stoc!(f, x, y, t)
        f[1] = +σ*x[1]
        return nothing
end
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), 0.0)

# Define ensemble quantities
Ne = 1000
N_esc = 0
N_uns = 0

# Define parameter values
Ny = 100
counter = 1
y_values = LinRange(0.0,3.0,Ny)

# Array for ensemble paths' endpoints 
Xt = fill(Float64[], 0) 
# Arrays for storing the ensemble's entire simple sample paths for a specific value of y
St = fill(Float64[], 0)
t = Float64[]

println("Simulating the ensemble sample paths")
# Loop over the ensemble sample paths
@showprogress for y in y_values
        local xt = Float64[]
        for j in 1:Ne
                local normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), y)
                local sol = solve(normal_form, EM(), dt=δt)
                # Store the ensemble path endpoint 
                push!(xt, sol[1,end])
                # Store the y-valued entire sample path
                if counter == 34 
                        push!(St, sol[1,:])
                        if j == 1
                                global t = sol.t
                        end
                end
        end
        push!(Xt, xt)
        global counter += 1
end

# Post processing to get the ensemble endpoint matrix
Et = Array{Float64, 2}(undef, Ne, Ny)
for j in 1:Ny
        ensemble_at_y = Xt[j]
        for k in 1:Ne
                Et[k,j] = ensemble_at_y[k] 
        end
end

# Plot the ensemble's endpoints parametric paths
fig1 = Figure(resolution = (1200, 400))
ax1 = Axis(fig1[1,1:4],
    ylabel = L"x",
    limits = ((y_values[1],y_values[end]), nothing)
)
for j in 1:Ne
        lines!(ax1, y_values, Et[j,:], color = :grey, linewidth=1)
end

# Plot the ensemble's average parametric path
At = Vector{Float64}(undef, Ny)
for j in 1:Ny
        At[j] = mean(Et[:,j])
end
lines!(ax1, y_values, At, color = :blue, linewidth=2)

# Plot the ensemble's variance at specific values of y
Vt = Vector{Float64}(undef, Ny)
for j in 1:Ny
        Vt[j] = var(Et[:,j]; mean=At[j]) 
end
lines!(ax1, y_values, Vt, color = :red, linewidth=2)

# Plot the ensemble' endpoints distribution at specific y-values
ax3 = Axis(fig1[2,1],
    title = L"y=0.90"
)
hist!(ax3, Et[:,31], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)
ax4 = Axis(fig1[2,2],
    title = L"y=2.00"
)
hist!(ax4, Et[:,67], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)
ax5 = Axis(fig1[2,3],
    title = L"y=2.30"
)
hist!(ax5, Et[:,77], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)
ax6 = Axis(fig1[2,4],
    title = L"y=2.75"
)
hist!(ax6, Et[:,92], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

# Plot the (deterministic) critical manifold
y_pos = range(y_values[1], y_values[end], length=10)
stable = 0.0 .+ 1.0*y_pos
unstable = 0.0 .+ 0.0*y_pos
lines!(ax1, y_pos, stable, color = :black, linewidth = 1.5)
lines!(ax1, y_pos, unstable, color = :black, linewidth = 1.5, linestyle = :dash)

# Plot a single sample path of the ensemble for a single parameter value to check steady-state is reached
fig2 = Figure(resolution = (1200, 400))
ax2 = Axis(fig2[1,1],
    ylabel = L"x",
    title = L"y=1"
)
for j in 1:Ne
        lines!(ax2, t, St[j])
end

# Export the figures
save("./results/Ensemble.png", fig1)
save("./results/SamplePaths.png", fig2)

#=
# Define uninitialised array to store a single sample path 
xt = Vector{Float64}(undef, Nt)
# Define uninitialised matrix to store the ensemble entire time-series
Xt = Array{Float64, 2}(undef, Ne, Nt)
# Define uninitialised matrix to store the ensemble sample paths' final states
Et = Array{Float64, 2}(undef, Ny, Ne)
# Define uninitialised matrix to store the parametric, ensemble-averaged, sample paths
At = Array{Float64, 2}(undef, Ny, Nt)

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
ax = Axis(fig1[1,1:4],
    xlabel = L"\tau",
    ylabel = L"x",
    limits = ((0,T), nothing)
)
for n in 1:Ny
        lines!(ax, LinRange(0,T,Nt), At[n,:], linewidth = 1)
end

# Plot the parametrised ensemble and temporal averaged sample paths 
ax2 = Axis(fig1[2,1:4],
    xlabel = L"y",
    ylabel = L"x",
    xticks = [0,1,2,3],
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
y_pos = range(3,0,length=10)
stable_pos = 0.0 .+ 1.0*y_pos 
unstable_pos = 0.0 .+ 0.0*y_pos 
lines!(y_pos, stable_pos, color = :black, linewidth = 0.5)
lines!(y_pos, unstable_pos, color = :black, linewidth = 0.5, linestyle = :dash)
=#

#=
# Plot P-bifurcation point of the transcritical normal form
yP = (σ^2)/2.0
lines!(ax1, yP.*ones(Float16,10), range(0,3.5,length=10), color = :black, linewidth = 1, linestyle = :dash)
text!(ax1, L"y^{+}_{P}", position = (yP+0.05,2.5), align = (:center,:center), color = :black, fontsize = 20)
lines!(ax2, yP.*ones(Float16,10), range(0,3.5,length=10), color = :black, linewidth = 1, linestyle = :dash)
text!(ax2, L"y^{+}_{P}", position = (yP+0.05,2.5), align = (:center,:center), color = :black, fontsize = 20)
=#

# Plot distributions of parametrised ensemble sample paths at the last state 
#=
ax3 = Axis(fig1[3,1],
         title = L"y=y^{+}_P",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax3, Et[11,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)
ax31 = Axis(fig1[4,1],
         title = L"y=y^{+}_P",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax31, At[11,:], bins = 50, color = :blue, strokecolor = :black, strokewidth = 1)

ax4 = Axis(fig1[3,2],
         title = L"y=1",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax4, Et[32,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)
ax41 = Axis(fig1[4,2],
         title = L"y=y^{+}_P",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax41, At[32,:], bins = 50, color = :blue, strokecolor = :black, strokewidth = 1)

ax5 = Axis(fig1[3,3],
         title = L"y=2",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax5, Et[66,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)
ax51 = Axis(fig1[4,3],
         title = L"y=y^{+}_P",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax51, At[66,:], bins = 50, color = :blue, strokecolor = :black, strokewidth = 1)

ax6 = Axis(fig1[3,4],
         title = L"y\,\lesssim\,3",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax6, Et[99,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)
ax61 = Axis(fig1[4,4],
         title = L"y=y^{+}_P",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}",
         limits = (0,nothing,nothing,nothing)
)
hist!(ax61, At[99,:], bins = 50, color = :blue, strokecolor = :black, strokewidth = 1)

# Plot the theoretical and numerical variance of the parametrised distribution 
fig2 = Figure(resolution = (800, 600))
ax1 = Axis(fig2[1,1],
           xlabel = L"y",
           ylabel = L"v_s",
           limits = ((0.0,3.0), (0,8))
)
# Numerical values computed from time-series
VarE = zeros(Float64, Ny)
VarA = zeros(Float64, Ny)
VarD = zeros(Float64, Ny)
for n in 1:Ny
        VarE[n] = varm(Et[n,:], mean(Et[n,:]); corrected=false)
        VarA[n] = mean((Et[n,:]).^2.0) - (mean(Et[n,:]))^2.0
        term = 0.0
        for m in 1:Ne
                term = term + (Et[n,m] - mean(Et[n,:]))^2.0
        end
        VarD[n] = (1/Ne)*term 
end
l1 = lines!(ax1, y_values, VarE, color = :red)
l2 = lines!(ax1, y_values, VarA, color = :blue, linestyle = :dash)
l3 = lines!(ax1, y_values, VarD, color = :green, linestyle = :dashdot)
# Position of the P-bifurcation
lines!(yP.*ones(Float16,100), range(0,8,length=100), color = :black, linewidth = 1.5, linestyle = :dash)
text!(L"y^{+}_{P}", position = (yP+0.1,5.0), align = (:center,:center), color = :black, fontsize = 20)
# Theoretical values computed from the steady-state Fokker-Planck distribution in closed form [Kuehn, 2011] 
VarK = zeros(Float64, Ny)
for n in 1:Ny
        VarK[n] = (y_values[n])^2.0 + 
                  (y_values[n]*σ^2)/2.0 -
                  (4.0^(-y_values[n]/σ^2.0))*y_values[n]*((1/σ^2.0)^((-2.0*y_values[n])/σ^2.0))*(σ^2.0)*gamma(1.0+(2*y_values[n])/σ^2.0) +
                  (2.0^((-4.0*y_values[n])/σ^2.0))*((y_values[n])^2.0)*((1/σ^2.0)^((-4.0*y_values[n])/σ^2.0))*(gamma((2.0*y_values[n])/σ^2.0))^(2.0)
end
l4 = lines!(ax1, y_values, VarK, color = :purple)
=#

# Export the figures
#save("SteadyStateFokkerPlanck.png", fig1)
#save("Variance.png", fig2)
