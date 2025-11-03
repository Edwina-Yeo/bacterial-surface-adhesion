#=
main:
- Julia version: 
- Author: edwinayeo
- Date: 2025-01-09
=#
using Random
using LinearAlgebra
using Printf
using DelimitedFiles
using Trapz
using Dierckx
using Plots
using StatsBase
using Distributions

# using BenchmarkTools


# function which takes in sol=(x,y,theta) of all current particles and returns, xb,yb,tb, if they have hit the bottom and sol=(x,y,theta) of those whcih havent
function swim(t, Vs, beta, Pe_r,Lx, Ly, dt)

        # step 1: update all particle positions
		master_sol_x[active_particles]=master_sol_x[active_particles] .+ (dt * (master_sol_y[active_particles] .+  Vs * cos.(master_sol_theta[active_particles])))
        master_sol_y[active_particles]= master_sol_y[active_particles] .+ (dt * (Vs * sin.(master_sol_theta[active_particles])))
        master_sol_theta[active_particles] =master_sol_theta[active_particles] .+ dt * ( - 1 / 2 .+ beta / 2 * cos.(2 * master_sol_theta[active_particles]) ).+  sqrt(2 * dt / Pe_r) * (randn(length(master_sol_theta[active_particles])))

        # now apply the no flux BC on top:
        master_sol_y[active_particles] = [yi > Ly ? Ly-(yi-Ly) : yi for yi in master_sol_y[active_particles]]


		# save the locations at which they hit the bottom
		master_sol_xb[active_particles .& (master_sol_y .< 0)] .= master_sol_x[active_particles .& (master_sol_y .<= 0)]
		master_sol_tb[active_particles .& (master_sol_y .< 0)] .= t

		# set indices of those which hit the bottom to be zero.
        # print(master_sol_y[active_particles] .<= 0)
		active_particles[active_particles] .= master_sol_y[active_particles] .> 0

        #set indices of those which exit on right to be zero.
		active_particles[active_particles] .= master_sol_x[active_particles] .<= Lx



    return
end

# diffusive jeffrey orbit: must be generated using series solution to PDE.
function calc_expansion_sol(Per, beta, theta_vec)
    N = 100 # number of terms from series solution
    a0 = 1 / (2 * π)

    A = zeros(2N, 2N)
    B = zeros(2N)

    # First equation
    A[1, 1] = -2
    A[1, 2] = beta
    A[1, N+1] = -8 / Per
    B[1] = -2 * beta * a0

    # Equations for a_k
    for i in 2:N-1
        A[i, i] = -2
        A[i, i-1] = beta
        A[i, i+1] = beta
        A[i, i+N] = -8 / Per * i
    end

    A[N, N] = -2
    A[N, N-1] = beta
    A[N, 2N] = -8 / Per * N

    # Equations for b_k
    A[N+1, N+1] = 2
    A[N+1, N+2] = -beta
    A[N+1, 1] = -8 / Per

    for i in N+2:2N-1
        A[i, i] = 2
        A[i, i-1] = -beta
        A[i, i+1] = -beta
        A[i, i-N] = -8 / Per * (i - N)
    end

    A[2N, 2N] = 2
    A[2N, 2N-1] = -beta
    A[2N, N] = -8 / Per * N

    coef_vec = A \ B

    sol = a0 * ones(length(theta_vec))
    for i in 1:N
        sol .+= coef_vec[i] * cos.(2 * i * theta_vec)
    end
    for i in N+1:2N
        sol .+= coef_vec[i] * sin.(2 * (i - N) * theta_vec)
    end

    return sol
end

 beta = parse(Float64, ARGS[1])
 Vs=parse(Float64,ARGS[2])
 Pe_r =parse(Float64,ARGS[3])

#beta=0.1
#Vs=1.68654808542314
#Pe_r=00790569415042095



Np = 2
dt = 0.01
T = 50
N_timesteps=Int(T/dt)
t = 0
Ly = 1.5
Lx = 3

Nbx = 150
Nby = 200
Nt = 50
L_tol = 0.1
dt_save=10
global N_arrive = 100





size_array=Int(N_arrive*N_timesteps)# size of all the space we need
global master_sol_x=zeros(size_array,) # for solution x

global master_sol_xb=zeros(size_array,) # for stuck bacteria
global master_sol_tb=zeros(size_array,) # time at which they stuck

# preallocate all random initial positions
global master_sol_y=sqrt.(rand(size_array,)).*Ly

# preallocate all random initial thetas
if beta>0
    # generate universal spline function to be called by solver
    theta_inits= range(0, 2π, length=500)
    solution = calc_expansion_sol(Pe_r, beta, theta_inits)
    cdf = cumsum(solution)
    jeffrey_spline = Spline1D(cdf, theta_inits; k=3)
    global master_sol_theta=jeffrey_spline(rand(size_array,))
    print("here",size(master_sol_x),size(master_sol_theta))
else
    global master_sol_theta=2*pi.*rand(size_array,)

end

# create logical array of active particles
global active_particles = falses(size_array,) # Initialize all particles as active



# Function to add new particles at the inlet
function add_particles!(active_particles,last_index,N_arrive)
    global active_particles[(last_index-1)*(N_arrive)+1:last_index*N_arrive] .= true # Mark new particles as active

return 
end

@time begin

# Start timing
global i_step=1
global t=0

while i_step<N_timesteps

swim(t,Vs, beta, Pe_r, Lx, Ly, dt) # progress particles forwards
add_particles!(active_particles,i_step,N_arrive) # add particles

 global i_step+=1
global t+=dt

end
# master_sol_xb[master_sol_tb.>0]
writedlm("data/Vs$Vs-beta-$beta-Per-$Pe_r-xbs.txt", master_sol_xb[master_sol_tb.>0], "   ")
writedlm("data/Vs$Vs-beta-$beta-Per-$Pe_r-tbs.txt", master_sol_tb[master_sol_tb.>0], "   ")


end