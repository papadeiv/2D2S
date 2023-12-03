using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Colors
using Statistics, ProgressMeter 

CairoMakie.activate!()
printstyled("PITCHFORK NORMAL FORM\n"; bold=true, underline=true, color=:light_blue)

# Define the parameters of the process
σ = 0.1
ε = 0.02

# Define initial states
x0 = [0.0,-1.0]

# Define temporal evolution quantities
T = 50.00
δt = 1e-3

# Define the transcritical form
function iip_det!(f, x, y, t)
        f[1] = x[2]*x[1] + (x[1])^3
        f[2] = ε
        return nothing
end
function iip_stoc!(f, x, y, t)
        f[1] = +σ
        f[2] = 0.0
        return nothing
end
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))

# Define ensemble quantities
Ne = 1000
N_esc = 0
N_uns = 0
critical_y = Float64[]

# Define parameter values
Ny = 100
y_values = LinRange(x0[2],0,Ny)

# Define the 'trajectory is escaped' termination condition
function condition_1(x, t, integrator)
        x[1] < -sqrt(abs(x[2])) || x[1] > +sqrt(abs(x[2]))
end
function condition_2(x, t, integrator)
        x[1] < 5*x[2] - 0.5
end
function condition_3(x, t, integrator)
        x[1] < -1.0 
end
function affect!(integrator)
        terminate!(integrator)
end
cb_1 = DiscreteCallback(condition_1, affect!, save_positions=(false,false))
cb_2 = DiscreteCallback(condition_2, affect!, save_positions=(false,false))
cb_3 = DiscreteCallback(condition_3, affect!, save_positions=(false,false))

# Define arrays to store the ensemble sample paths for escaped and unescaped trajectories 
xt = fill(Float64[], 0) 
yt = fill(Float64[], 0)
xt_esc = fill(Float64[], 0)
yt_esc = fill(Float64[], 0)
max_yt_series = Float64[]

println("Simulating the ensemble sample paths")
# Loop over the ensemble sample paths
@showprogress for j in 1:Ne
        local sol = solve(normal_form, EM(), dt=δt, callback=cb_1, verbose=false)
        if sol[1,end] < -sqrt(abs(sol[2,end])) || sol[1,end] > +sqrt(abs(sol[2,end])) # It means that this trajectory has escaped the manifold 
                # Store the times series up to the escape
                push!(xt_esc, sol[1,:])
                push!(yt_esc, sol[2,:])
                # Store the value of the slow variable at the unstable critical submanifold
                push!(critical_y, sol[2,end])
                # Update the counter of escaped trajectories
                global N_esc = N_esc + 1
                # Check if the sample path has the longest time series so far and in that case update the corresponding variable
                if size(sol[2,:],1) > size(max_yt_series,1)
                        empty!(max_yt_series)
                        global max_yt_series = sol[2,:]
                end
        else # It means that the trajectory stayed bounded for the whole time
                push!(xt, sol[1,:])
                push!(yt, sol[2,:])
                # Update the counter of unescaped trajectories
                global N_uns = N_uns + 1
        end
end
# Plot the ensemble's sample paths
fig1 = Figure(resolution = (1200, 400))
ax1 = Axis(fig1[1,1:2],
    ylabel = L"x",
    limits = ((x0[2],y_values[end]), (-1,+1))
)
ax2 = Axis(fig1[2,1:2],
    ylabel = L"x",
    limits = ((x0[2],y_values[end]), (-1,+1))
)
# Unescaped trajectories
for n in 1:N_uns # == size(xt,1)
        lines!(ax1, yt[n], xt[n], linewidth = 0.5, color = :gray)
end
# Escaped trajectories
for n in 1:N_esc # == size(xt_esc,1)
        lines!(ax2, yt_esc[n], xt_esc[n], linewidth = 1, color = :red)
end
# Plot the (deterministic) critical manifold
stable = 0.0.*y_values 
unstable1 = sqrt.(-y_values) 
unstable2 = -unstable1
lines!(ax1, y_values, stable, color = :black, linewidth = 1.5)
lines!(ax1, y_values, unstable1, color = :black, linewidth = 1.5, linestyle = :dash)
lines!(ax1, y_values, unstable2, color = :black, linewidth = 1.5, linestyle = :dash)
# Plot percentage of escaped trajectories at different y values
ax3 = Axis(fig1[1:2,3],
         xlabel = L"y",
         ylabel = "Escaped traj.",
         limits = ((-1,0), nothing)
)
pct_esc = Float64[]
push!(pct_esc, 0.0)
counter = 0
println("Computing the parameter distribution of the escaped trajectories")
# Loop over the discretised parameter values
@showprogress for n in 2:Ny
        for m in 1:N_esc
                w = yt_esc[m]
                if w[end] > y_values[n-1] && w[end] < y_values[n] 
                        global counter = counter + 1
                end
        end
        push!(pct_esc, (counter/Ne)*100)
end
scatter!(ax3, y_values, pct_esc, color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 7)
text!(-0.75, 80, text="Tot. esc. = " .*string.(round((N_esc/Ne)*100, digits=2)) .* "%", align = (:center, :center))
# Plot the ensemble moments of the distributions for different parameter values
M = Float64[]
V = Float64[]
println("Computing the variance of the (unescaped) ensemble paths")
# Loop over the parameter values (timesteps)
@showprogress for n in 1:size(max_yt_series,1) 
        # Compute the ensemble mean
        MEAN = 0.0
        # Loop over the escaped ensemble paths all at the n-th timestep
        N_esc_crit = 0
        for m in 1:N_esc # == size(xt_esc,1)
                u = xt_esc[m]
                w = yt_esc[m]
                if n <= size(w,1)
                        if w[n] < critical_y[m] 
                                MEAN = MEAN + u[n]
                                N_esc_crit = N_esc_crit + 1
                        end
                end
        end
        MEAN = MEAN/(N_uns + N_esc_crit)
        # Compute the ensemble variance
        VAR = 0.0
        # Loop over the escaped ensemble paths all at the n-th timestep
        N_esc_crit = 0
        for m in 1:N_esc # == size(xt_esc,1)
                u = xt_esc[m]
                w = yt_esc[m]
                if n <= size(w,1)
                        if w[n] < critical_y[m] 
                                VAR = VAR + (u[n] - MEAN)^(2.0)
                                N_esc_crit = N_esc_crit + 1
                        end
                end
        end
        VAR = VAR/(N_uns + N_esc_crit)
        # Store the moments of the ensemble at the n-th timestep
        push!(M, MEAN)
        push!(V, VAR)
end
lines!(ax1, max_yt_series, M, color = :blue, linewidth = 1.5)
fig2 = Figure(resolution = (1200, 1200))
ax4 = Axis(fig2[1,1],
    xlabel = L"y",
    ylabel = "Variance",
    limits = ((-1,0), nothing)
)
lines!(ax4, max_yt_series, V, color = :black, linewidth = 1.5)
# Plot the theoretical variance of the steady-state Fokker-Planck distribution of the normal form
theoretical_y = LinRange(x0[2],0,200)
theoretical_V = [0.0050782,0.00510455,0.00513117,0.00515808,0.00518528,0.00521277,0.00524056,0.00526865,0.00529706,0.00532577,0.0053548,0.00538416,0.00541385,0.00544388,0.00547424,0.00550495,0.00553602,0.00556745,0.00559924,0.0056314,0.00566395,0.00569688,0.00573021,0.00576393,0.00579807,0.00583262,0.00586759,0.005903,0.00593885,0.00597514,0.0060119,0.00604911,0.00608681,0.00612499,0.00616366,0.00620284,0.00624253,0.00628275,0.0063235,0.0063648,0.00640666,0.00644909,0.0064921,0.00653571,0.00657992,0.00662476,0.00667022,0.00671634,0.00676312,0.00681058,0.00685873,0.00690759,0.00695717,0.0070075,0.00705859,0.00711045,0.00716311,0.00721659,0.0072709,0.00732606,0.00738211,0.00743905,0.00749691,0.00755572,0.00761551,0.00767628,0.00773808,0.00780093,0.00786485,0.00792988,0.00799605,0.00806339,0.00813193,0.0082017,0.00827275,0.00834511,0.00841882,0.00849392,0.00857045,0.00864845,0.00872798,0.00880908,0.0088918,0.0089762,0.00906233,0.00915025,0.00924001,0.0093317,0.00942537,0.00952109,0.00961895,0.00971901,0.00982138,0.00992613,0.0100334,0.0101432,0.0102557,0.010371,0.0104892,0.0106105,0.010735,0.0108629,0.0109942,0.0111292,0.0112681,0.011411,0.0115582,0.0117098,0.0118662,0.0120276,0.0121942,0.0123664,0.0125445,0.0127288,0.0129197,0.0131175,0.0133228,0.0135358,0.0137572,0.0139874,0.0142268,0.0144761,0.0147359,0.0150066,0.0152889,0.0155834,0.0158908,0.0162117,0.0165466,0.0168962,0.0172611,0.0176419,0.018039,0.0184529,0.018884,0.0193325,0.0197986,0.0202824,0.0207838,0.0213024,0.0218379,0.0223895,0.0229565,0.0235377,0.0241317,0.0247369,0.0253514,0.025973,0.0265994,0.0272276,0.0278549,0.0284779,0.0290931,0.0296968,0.0302852,0.0308542,0.0313998,0.0319175,0.0324033,0.0328528,0.0332617,0.0336261,0.033942,0.0342054,0.0344128,0.034561,0.0346468,0.0346675,0.0346208,0.0345044,0.0343168,0.0340567,0.033723,0.0333152,0.0328331,0.0322768,0.0316468,0.0309439,0.0301693,0.0293244,0.0284109,0.0274308,0.0263862,0.0252796,0.0241136,0.022891,0.0216147,0.0202877,0.0189134,0.017495,0.0160359,0.0145395,0.0130096,0.0114496,0.00986316,0.00825409,0.00662607,0.00498287,0.00332826,0.00166603]
lines!(ax4, theoretical_y, theoretical_V, color = :blue, linewidth = 3.5)
# Export figures
save("../results/StocPitchBifDistributions.png", fig1)
save("../results/StocPitchBifVariance.png", fig2)
