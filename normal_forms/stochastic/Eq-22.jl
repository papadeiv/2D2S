using DifferentialEquations
using LaTeXStrings, CairoMakie
using Statistics


#     Ornstein-Uhlenbeck process
# (linear drift, constant diffusion)


# Define the parameters of the process
σ = 0.1 
ε = 0.02
α = 1.0
trueVar = 0.005
∂N = 0.35

# Define the deterministic dynamics
function iip_det!(f, x, y, t)
        f[1] = -(α/ε)*x[1] 
        f[2] = 1.0 
        return nothing
end
# Define the stochastic dynamics
function iip_stoc!(f, x, y, t)
        f[1] = + σ/sqrt(ε)
        f[2] = 0.0
        return nothing
end

# Define the evolution parameters 
x0 = [0.0, -2.0]
T = 2.00
δt = 1e-2

# Build the (non-parametric) dynamical system
OUP = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))

# Evolve a sample path of the stochastic process 
Xt = solve(OUP, EM(), dt=δt)

# Print the sample path 
#=
for n in 1:size(Xt,2)
        printstyled("Solution x=", Xt[n], " at time t=", Xt.t[n], "\n"; bold=true, underline=true, color=:light_green)
end
=#

# Plot the sample path 
CairoMakie.activate!()
fig1 = Figure(resolution = (1000, 400))
ax = Axis(fig1[1,2],
    xlabel = L"\tau",
    ylabel = L"x",
    yticks = [-0.5,-0.35,0,0.35,0.5] 
)
limits!(ax,x0[2],x0[2]+T,-∂N-0.15,∂N+0.15)
lines!(ax, Xt[2,:], Xt[1,:], color= :black, linewidth = 1)
y_var = range(x0[2],x0[2]+T,length=1000) 
x_var = trueVar .+ 0.0*y_var 
x_∂N = ∂N .+ 0.0*y_var
lines!(y_var, x_var, color = :red, linewidth = 1, linestyle = :dash)
lines!(y_var, -x_var, color = :red, linewidth = 1, linestyle = :dash)
lines!(y_var, x_∂N, color = :blue, linewidth = 1, linestyle = :dash)
lines!(y_var, -x_∂N, color = :blue, linewidth = 1, linestyle = :dash)

# Plot distribution of the time-series
ax2 = Axis(fig1[1,1],
         aspect = 0.75,
         xlabel = L"x_{\tau}",
         ylabel = L"n^{°} realizations",
         xticks = [0,0.04],
         yticks = [0]
)
hist!(ax2, Xt[1,:], bins = 50, normalization = :probability, color = :red, strokecolor = :black, strokewidth = 1, direction = :x)

# Plot the time-varying moments of the distribution 
ax3 = Axis(fig1[2,2],
           xlabel = L"\tau",
)
time = LinRange(0,T,size(Xt,2)) 
# Numerical values computed from time-series
E = mean(Xt[1,:]) .+ 0.0*time
Var = zeros(Float64, size(Xt,2))
for n in 1:size(Xt,2)
        Var[n] = var(Xt[1,1:n]; mean = E[1])
end
#println(last(Var))
# Theoretical values [Kuehn, 2011]
E_exact = x0[1].*exp.((-α.*time)./ε)
Var_exact = (x0[1] - (σ^2)/(2*α)).*exp.((-2*α.*time)./ε) .+ (σ^2)/(2*α)
#println(last(Var_exact))
limits!(ax3,0,T,0,0.01)
l1 = lines!(ax3, time, Var, color = :red)
l2 = lines!(ax3, time, Var_exact, color = :blue, linestyle = :dash)
Legend(fig1[2,1], [l1,l2], [L"Var_{num}", L"Var_{exact}"])

# Export the figure
save("OUP.png", fig1)

# Create animation
points = Observable(Point2f[(x0[2],x0[1])])
fig2 = Figure(resolution = (500, 200))
ax3 = Axis(fig2[1,1],
           xticks = [1], 
          yticks = [1] 
)
limits!(ax,-2,0,-0.5,0.5)
lines!(points, color = :black)
y_var = range(-2,0,length=1000) 
x_var = 0.01 .+ 0.0*y_var 
x_deltaN = ∂N .+ 0.0*y_var
lines!(y_var, x_var, color = :red, linewidth = 1, linestyle = :dash)
lines!(y_var, -x_var, color = :red, linewidth = 1, linestyle = :dash)
lines!(y_var, x_deltaN, color = :blue, linewidth = 1, linestyle = :dash)
lines!(y_var, -x_deltaN, color = :blue, linewidth = 1, linestyle = :dash)

frames = 1:size(Xt,2)
framerate = 40

record(fig2, "Equation 22.gif", frames;
       framerate=framerate) do frame 
        println(frame)
       new_point = Point2f(Xt[2,frame],Xt[1,frame])
       points[] = push!(points[], new_point)
end
