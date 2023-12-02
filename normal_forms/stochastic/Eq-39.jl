using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Colors
using Statistics, SpecialFunctions


# Normal forms stochastic bifurctions 

# Define the parameters of the process
σ = 0.1

# Define the fold normal form 
function iip_fold_det!(f, x, y, t)
        f[1] = -y[1] - (x[1])^2
        return nothing
end
function iip_fold_stoc!(f, x, y, t)
        f[1] = +σ
        return nothing
end

# Define the transcritical form
function iip_trans_det!(f, x, y, t)
        f[1] = y[1]*x[1] - (x[1])^2
        return nothing
end
function iip_trans_stoc!(f, x, y, t)
        f[1] = +σ
        return nothing
end

# Define the pitchfork normal form 
function iip_pitch_det!(f, x, y, t)
        f[1] = y[1]*x[1] + (x[1])^3
        return nothing
end
function iip_pitch_stoc!(f, x, y, t)
        f[1] = +σ
        return nothing
end

# Define initial states
x0_SP = [1.5]
x0_TC = [0.1]
x0_PF = [0.5]

# Define ensemble quantities
Ne = 1000

# Define parameter values
Ny = 100
counter = 1
y_values = LinRange(-1,0,Ny)

# Define temporal evolution quantities
T_SP = 4.00
T_TC = 3.00
T_PF = 1.00
δt = 1e-3
Nt_SP = ceil(Int64, T_SP/δt) + 2
Nt_TC = ceil(Int64, T_TC/δt) + 2
Nt_PF = ceil(Int64, T_PF/δt) + 1

# Define uninitialised array to store a single sample path 
xt_SP = Vector{Float64}(undef, Nt_SP)
xt_TC = Vector{Float64}(undef, Nt_TC)
xt_PF = Vector{Float64}(undef, Nt_PF)
# Define uninitialised matrix to store the ensemble entire time-series
Xt_SP = Array{Float64, 2}(undef, Ne, Nt_SP)
Xt_TC = Array{Float64, 2}(undef, Ne, Nt_TC)
Xt_PF = Array{Float64, 2}(undef, Ne, Nt_PF)
# Define uninitialised matrix to store the ensemble sample paths' final states
Et_SP = Array{Float64, 2}(undef, Ny, Ne)
Et_TC = Array{Float64, 2}(undef, Ny, Ne)
Et_PF = Array{Float64, 2}(undef, Ny, Ne)
# Define uninitialised matrix to store the parametric, ensemble-averaged, sample paths
At_SP = Array{Float64, 2}(undef, Ny, Nt_SP)
At_TC = Array{Float64, 2}(undef, Ny, Nt_TC)
At_PF = Array{Float64, 2}(undef, Ny, Nt_PF)

# FOLD BIFURCATION 
println("Fold bifurcation")
for y in y_values
        printstyled(counter, ") Ensemble simulation for y=", y, "\n"; bold=true, underline=true, color=:light_green)
        # Perform an ensemble simulation of Ne sample paths for the current parameter value y
        for j in 1:Ne
                # Build the (non-parametric) dynamical system
                local AandB = SDEProblem(iip_fold_det!, iip_fold_stoc!, x0_SP, (0.0, T_SP), y)
                # Evolve a sample path of the stochastic process for the current parameter value 
                sol = solve(AandB, EM(), dt=δt)
                # Store the sample path into the appropriate trajectory matrix 
                global Xt_SP[j,:] = sol[1,:]
        end
        # Compute the ensemble averaged trajectory for the current parameter value y
        for k in 1:Nt_SP
                global xt_SP[k] = mean(Xt_SP[:,k]) 
        end
        # Store the ensemble's sample paths' final states into the appropriate trajectory matrix
        global Et_SP[counter,:] = Xt_SP[:,end]
        # Store the sample path into the appropriate trajectory matrix
        global At_SP[counter,:] = xt_SP 
        printstyled("   Ensemble avg=", mean(Et_SP[counter,:]), "\n"; color=:light_blue)
        printstyled("   Temporal avg=", mean(At_SP[counter,:]), "\n"; color=:light_red)
        # Update the counter
        global counter += 1
end

# Plot the parametrised, ensemble-averaged, sample paths 
CairoMakie.activate!()
fig1 = Figure(resolution = (1200, 900))
ax = Axis(fig1[1,1:3],
    xlabel = L"\tau",
    ylabel = L"x",
    limits = ((0,T_SP), (0,1.55))
)
for n in 1:Ny
        lines!(ax, LinRange(0,T_SP,Nt_SP), At_SP[n,:], linewidth = 1)
end

# Plot the parametrised ensemble and temporal averaged sample paths 
ax2 = Axis(fig1[2,1:3],
    xlabel = L"y",
    ylabel = L"x",
    xticks = [-0.75,-0.5,-0.25],
    limits = ((y_values[1],last(y_values)), nothing)
)
# The ensemble average is computed over the final states of the sample paths of the ensemble
nsmb_avg = zeros(Ny)
# The temporal average is computed over the ensemble's averaged states across the time series
temp_avg = zeros(Ny)
for n in 1:Ny
        nsmb_avg[n] = mean(Et_SP[n,:])
        temp_avg[n] = mean(At_SP[n,:])
end
lines!(ax2, y_values, nsmb_avg, color = :red, linewidth = 1)
lines!(ax2, y_values, temp_avg, color = :blue, linewidth = 1)

# Plot the (deterministic) critical manifold
stable = sqrt.(-y_values)
unstable = -sqrt.(-y_values)
lines!(y_values, stable, color = :black, linewidth = 0.5)
lines!(y_values, unstable, color = :black, linewidth = 0.5, linestyle = :dash)

# Plot distributions of parametrised ensemble sample paths at the last state 
ax3 = Axis(fig1[3,1],
         title = L"y=-0.75",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}"
)
hist!(ax3, Et_SP[26,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

ax4 = Axis(fig1[3,2],
         title = L"y=-0.5"
)
hist!(ax4, Et_SP[51,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

ax5 = Axis(fig1[3,3],
         title = L"y=-0.25"
)
hist!(ax5, Et_SP[75,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

# Export the figures
save("Fold.png", fig1)

# TRANSCRITICAL BIFURCATION 
counter = 1
println("Transcritical bifurcation")
for y in y_values
        printstyled(counter, ") Ensemble simulation for y=", y, "\n"; bold=true, underline=true, color=:light_green)
        # Perform an ensemble simulation of Ne sample paths for the current parameter value y
        for j in 1:Ne
                # Build the (non-parametric) dynamical system
                local AandB = SDEProblem(iip_trans_det!, iip_trans_stoc!, x0_TC, (0.0, T_TC), y)
                # Evolve a sample path of the stochastic process for the current parameter value 
                sol = solve(AandB, EM(), dt=δt)
                # Store the sample path into the appropriate trajectory matrix 
                global Xt_TC[j,:] = sol[1,:]
        end
        # Compute the ensemble averaged trajectory for the current parameter value y
        for k in 1:Nt_TC
                global xt_TC[k] = mean(Xt_TC[:,k]) 
        end
        # Store the ensemble's sample paths' final states into the appropriate trajectory matrix
        global Et_TC[counter,:] = Xt_TC[:,end]
        # Store the sample path into the appropriate trajectory matrix
        global At_TC[counter,:] = xt_TC 
        printstyled("   Ensemble avg=", mean(Et_TC[counter,:]), "\n"; color=:light_blue)
        printstyled("   Temporal avg=", mean(At_TC[counter,:]), "\n"; color=:light_red)
        # Update the counter
        global counter += 1
end

# Plot the parametrised, ensemble-averaged, sample paths 
CairoMakie.activate!()
fig2 = Figure(resolution = (1200, 900))
ax6 = Axis(fig2[1,1:3],
    xlabel = L"\tau",
    ylabel = L"x",
    limits = ((0,T_TC), nothing)
)
for n in 1:Ny
        lines!(ax6, LinRange(0,T_TC,Nt_TC), At_TC[n,:], linewidth = 1)
end

# Plot the parametrised ensemble and temporal averaged sample paths 
ax7 = Axis(fig2[2,1:3],
    xlabel = L"y",
    ylabel = L"x",
    xticks = [-0.75,-0.5,-0.25],
    limits = ((y_values[1],last(y_values)), nothing)
)
# The ensemble average is computed over the final states of the sample paths of the ensemble
nsmb_avg = zeros(Ny)
# The temporal average is computed over the ensemble's averaged states across the time series
temp_avg = zeros(Ny)
for n in 1:Ny
        nsmb_avg[n] = mean(Et_TC[n,:])
        temp_avg[n] = mean(At_TC[n,:])
end
lines!(ax7, y_values, nsmb_avg, color = :red, linewidth = 1)
lines!(ax7, y_values, temp_avg, color = :blue, linewidth = 1)

# Plot the (deterministic) critical manifold
y_neg = range(-1,0,length=10)
stable_neg = 0.0 .+ 0.0*y_neg
unstable_neg = 0.0 .+ 1.0*y_neg
lines!(y_neg, stable_neg, color = :black, linewidth = 0.5)
lines!(y_neg, unstable_neg, color = :black, linewidth = 0.5, linestyle = :dash)

# Plot distributions of parametrised ensemble sample paths at the last state 
ax8 = Axis(fig2[3,1],
         title = L"y=-0.75",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}"
)
hist!(ax8, Et_TC[26,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

ax9 = Axis(fig2[3,2],
         title = L"y=-0.5"
)
hist!(ax9, Et_TC[51,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

ax10 = Axis(fig2[3,3],
         title = L"y=-0.25"
)
hist!(ax10, Et_TC[75,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

# Export the figures
save("Transcritical.png", fig2)

# PITCHFORK BIFURCATION 
counter = 1
println("Pichfork bifurcation")
for y in y_values
        printstyled(counter, ") Ensemble simulation for y=", y, "\n"; bold=true, underline=true, color=:light_green)
        # Perform an ensemble simulation of Ne sample paths for the current parameter value y
        for j in 1:Ne
                # Build the (non-parametric) dynamical system
                local AandB = SDEProblem(iip_pitch_det!, iip_pitch_stoc!, x0_PF, (0.0, T_PF), y)
                # Evolve a sample path of the stochastic process for the current parameter value 
                sol = solve(AandB, EM(), dt=δt)
                # Store the sample path into the appropriate trajectory matrix 
                global Xt_PF[j,:] = sol[1,:]
        end
        # Compute the ensemble averaged trajectory for the current parameter value y
        for k in 1:Nt_PF
                global xt_PF[k] = mean(Xt_PF[:,k]) 
        end
        # Store the ensemble's sample paths' final states into the appropriate trajectory matrix
        global Et_PF[counter,:] = Xt_PF[:,end]
        # Store the sample path into the appropriate trajectory matrix
        global At_PF[counter,:] = xt_PF 
        printstyled("   Ensemble avg=", mean(Et_PF[counter,:]), "\n"; color=:light_blue)
        printstyled("   Temporal avg=", mean(At_PF[counter,:]), "\n"; color=:light_red)
        # Update the counter
        global counter += 1
end

# Plot the parametrised, ensemble-averaged, sample paths 
CairoMakie.activate!()
fig3 = Figure(resolution = (1200, 900))
ax11 = Axis(fig3[1,1:3],
    xlabel = L"\tau",
    ylabel = L"x",
    limits = ((0,T_PF), nothing)
)
for n in 1:Ny
        lines!(ax11, LinRange(0,T_PF,Nt_PF), At_PF[n,:], linewidth = 1)
end

# Plot the parametrised ensemble and temporal averaged sample paths 
ax12 = Axis(fig3[2,1:3],
    xlabel = L"y",
    ylabel = L"x",
    xticks = [-0.75,-0.5,-0.25],
    limits = ((y_values[1],last(y_values)), nothing)
)
# The ensemble average is computed over the final states of the sample paths of the ensemble
nsmb_avg = zeros(Ny)
# The temporal average is computed over the ensemble's averaged states across the time series
temp_avg = zeros(Ny)
for n in 1:Ny
        nsmb_avg[n] = mean(Et_PF[n,:])
        temp_avg[n] = mean(At_PF[n,:])
end
lines!(ax12, y_values, nsmb_avg, color = :red, linewidth = 1)
lines!(ax12, y_values, temp_avg, color = :blue, linewidth = 1)

# Plot the (deterministic) critical manifold
stable = 0.0.*y_values 
unstable1 = sqrt.(-y_values) 
unstable2 = -unstable1
lines!(y_values, stable, color = :black, linewidth = 0.5)
lines!(y_values, unstable1, color = :black, linewidth = 0.5, linestyle = :dash)
lines!(y_values, unstable2, color = :black, linewidth = 0.5, linestyle = :dash)

# Plot distributions of parametrised ensemble sample paths at the last state 
ax13 = Axis(fig3[3,1],
         title = L"y=-0.75",
         ylabel = L"p_s(x)",
         xlabel = L"x_{\tau}"
)
hist!(ax13, Et_PF[26,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

ax14 = Axis(fig3[3,2],
         title = L"y=-0.5"
)
hist!(ax14, Et_PF[51,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

ax15 = Axis(fig3[3,3],
         title = L"y=-0.25"
)
hist!(ax15, Et_PF[75,:], bins = 100, color = :red, strokecolor = :black, strokewidth = 1)

# Export the figures
save("Pitchfork.png", fig3)
