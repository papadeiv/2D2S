using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Colors
using Statistics, ProgressMeter 

CairoMakie.activate!()
printstyled("TRANSCRITICAL NORMAL FORM\n"; bold=true, underline=true, color=:light_blue)

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
        f[1] = x[2]*x[1] - (x[1])^2
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
        x[1] < x[2]
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

println("Simulating the ensemble sample paths")
# Loop over the ensemble sample paths
@showprogress for j in 1:Ne
        local sol = solve(normal_form, EM(), dt=δt, callback=cb_1, verbose=false)
        if sol[1,end] < sol[2,end] # It means that this trajectory has escaped the manifold 
                # Store the times series up to the escape
                push!(xt_esc, sol[1,:])
                push!(yt_esc, sol[2,:])
                # Store the value of the slow variable at the unstable critical submanifold
                push!(critical_y, sol[2,end])
                # Update the counter of escaped trajectories
                global N_esc = N_esc + 1
        else # It means that the trajectory stayed bounded for the whole time
                push!(xt, sol[1,:])
                push!(yt, sol[2,:])
                # Update the counter of unescaped trajectories
                global N_uns = N_uns + 1
        end
end
# Plot the ensemble's sample paths
fig1 = Figure(; size = (1200, 400), backgroundcolor = :transparent,)
ax1 = Axis(fig1[1,1:2],
    backgroundcolor = :transparent,
    ylabel = L"x",
    limits = ((x0[2],y_values[end]), (-1,+0.5))
)
ax2 = Axis(fig1[2,1:2],
    backgroundcolor = :transparent,
    ylabel = L"x",
    limits = ((x0[2],y_values[end]), (-1,+0.5))
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
y_neg = range(y_values[1],0,length=10)
stable_neg = 0.0 .+ 0.0*y_neg
unstable_neg = 0.0 .+ 1.0*y_neg
lines!(ax1, y_neg, stable_neg, color = :black, linewidth = 1.5)
lines!(ax1, y_neg, unstable_neg, color = :black, linewidth = 1.5, linestyle = :dash)
lines!(ax2, y_neg, stable_neg, color = :black, linewidth = 1.5)
lines!(ax2, y_neg, unstable_neg, color = :black, linewidth = 1.5, linestyle = :dash)
# Plot percentage of escaped trajectories at different y values
ax3 = Axis(fig1[1:2,3],
         backgroundcolor = :transparent,
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
@showprogress for n in 1:size(xt[end],1) 
        # Compute the ensemble mean
        MEAN = 0.0
        # Loop over the unescaped ensemble paths all at the n-th timestep
        for m in 1:N_uns # == size(xt,1)
                u = xt[m]
                MEAN = MEAN + u[n] 
        end
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
        # Loop over the unescaped ensemble paths all at the n-th timestep
        for m in 1:N_uns
                u = xt[m]
                VAR = VAR + (u[n] - MEAN)^(2.0) 
        end
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
lines!(ax1, yt[end], M, color = :blue, linewidth = 1.5)
fig2 = Figure(; size = (1200, 1200), backgroundcolor = :transparent)
ax4 = Axis(fig2[1,1],
    backgroundcolor = :transparent,
    xlabel = L"y",
    ylabel = "Variance",
    limits = ((-1,0), nothing)
)
lines!(ax4, yt[end], V, color = :black, linewidth = 3.5)
# Plot the theoretical variance of the steady-state Fokker-Planck distribution of the normal form
theoretical_y = LinRange(x0[2],0,200)
theoretical_V = [0.00510694,0.00513435,0.00516208,0.00519013,0.00521851,0.00524723,0.00527628,0.00530569,0.00533545,0.00536557,0.00539607,0.00542694,0.00545819,0.00548984,0.0055219,0.00555436,0.00558724,0.00562055,0.0056543,0.00568849,0.00572315,0.00575827,0.00579387,0.00582996,0.00586655,0.00590365,0.00594128,0.00597945,0.00601817,0.00605746,0.00609733,0.0061378,0.00617887,0.00622058,0.00626293,0.00630594,0.00634963,0.00639403,0.00643915,0.00648501,0.00653163,0.00657904,0.00662727,0.00667634,0.00672627,0.0067771,0.00682886,0.00688158,0.00693529,0.00699003,0.00704584,0.00710277,0.00716085,0.00722013,0.00728067,0.00734252,0.00740573,0.00747036,0.0075365,0.00760419,0.00767353,0.00774459,0.00781746,0.00789222,0.00796899,0.00804786,0.00812895,0.00821238,0.00829827,0.00838676,0.00847797,0.00857207,0.0086692,0.00876952,0.0088732,0.00898039,0.00909127,0.00920601,0.00932478,0.00944775,0.00957509,0.00970694,0.00984347,0.00998481,0.0101311,0.0102824,0.0104388,0.0106004,0.0107672,0.0109392,0.0111163,0.0112985,0.0114857,0.0116776,0.011874,0.0120747,0.0122793,0.0124875,0.0126988,0.0129127,0.0131286,0.0133459,0.0135641,0.0137823,0.0139999,0.0142161,0.01443,0.0146409,0.0148479,0.0150501,0.0152467,0.0154369,0.0156198,0.0157946,0.0159607,0.0161172,0.0162635,0.0163989,0.0165231,0.0166353,0.0167353,0.0168227,0.0168972,0.0169586,0.0170069,0.0170419,0.0170637,0.0170725,0.0170682,0.0170513,0.0170219,0.0169804,0.0169273,0.0168628,0.0167875,0.0167018,0.0166063,0.0165015,0.016388,0.0162662,0.0161368,0.0160004,0.0158575,0.0157086,0.0155543,0.0153951,0.0152316,0.0150642,0.0148935,0.0147199,0.0145438,0.0143656,0.0141859,0.0140049,0.0138229,0.0136405,0.0134577,0.013275,0.0130926,0.0129108,0.0127298,0.0125497,0.0123709,0.0121934,0.0120174,0.0118431,0.0116706,0.0115,0.0113314,0.0111649,0.0110005,0.0108384,0.0106786,0.0105211,0.010366,0.0102133,0.0100631,0.00991527,0.00976991,0.00962701,0.00948657,0.00934858,0.00921303,0.00907991,0.0089492,0.00882089,0.00869495,0.00857135,0.00845007,0.00833108,0.00821435,0.00809985,0.00798754,0.0078774,0.00776938,0.00766344,0.00755957,0.00745771,0.00735783,0.0072599]
lines!(ax4, theoretical_y, theoretical_V, color = :green, linewidth = 5)
# Export figures
save("../results/AddTransEnsemble.png", fig1)
save("../results/AddTransVariance.png", fig2)
