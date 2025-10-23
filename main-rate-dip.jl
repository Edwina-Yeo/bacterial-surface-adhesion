using Dierckx
const target = 1000000
#const target = 10000000 # t=1000


const block_size = 200
const dt = 0.01
const L_y = 1.5
const L_x = 3
cellsPerSecond = 10000


using Dierckx
const target = 10000000
const block_size = 200
const dt = 0.01
#const beta = 0
const L_y = 1.5
const L_x = 3
cellsPerSecond = 10000

using DelimitedFiles

print("effective time is ",target/(cellsPerSecond) )



gammas = exp10.(range(-1.5,2,20))

beta = parse(Float64, ARGS[1])
Vs_val=parse(Float64,ARGS[2])
D_r =parse(Float64,ARGS[3])
kappa =parse(Float64,ARGS[4])
rep =parse(Float64,ARGS[5])

#println(beta,Vs_val,D_r,kappa,rep)

 
 

L=750 # lengthscale for re-dimensionalising 
U=gammas*L # flow speeds
Vss=Vs_val./U 
Pers=U./(D_r*L);
alpha=31.83

dips=alpha/(L^2*U)


function run_sim(Pe_r, Vs, gamma, kappa,dip)

    D_step = sqrt(2 * dt / Pe_r)

    swim_x = zeros(block_size, 1)
    swim_y = zeros(block_size, 1)
    swim_theta = zeros(block_size, 1)
    swim_tstart = zeros(block_size, 1)

    stick_events = Tuple{Float64, Float64}[]

    # Prepare a spline for sampling from the Jeffrey distribution
    # diffusive jeffrey orbit: must be generated using series solution to PDE.
    function calc_expansion_sol(theta)
        N = 100 # number of terms from series solution
        a0 = 1 / (2 * π)

        A = zeros(2N, 2N)
        B = zeros(2N)

        # First equation
        A[1, 1] = -2
        A[1, 2] = beta
        A[1, N+1] = -8 / Pe_r
        B[1] = -2 * beta * a0

        # Equations for a_k
        for i in 2:N-1
            A[i, i] = -2
            A[i, i-1] = beta
            A[i, i+1] = beta
            A[i, i+N] = -8 / Pe_r * i
        end

        A[N, N] = -2
        A[N, N-1] = beta
        A[N, 2N] = -8 / Pe_r * N

        # Equations for b_k
        A[N+1, N+1] = 2
        A[N+1, N+2] = -beta
        A[N+1, 1] = -8 / Pe_r

        for i in N+2:2N-1
            A[i, i] = 2
            A[i, i-1] = -beta
            A[i, i+1] = -beta
            A[i, i-N] = -8 / Pe_r * (i - N)
        end

        A[2N, 2N] = 2
        A[2N, 2N-1] = -beta
        A[2N, N] = -8 / Pe_r * N

        coef_vec = A \ B

        sol = a0 * ones(length(theta))
        for i in 1:N
            sol .+= coef_vec[i] * cos.(2 * i * theta)
        end
        for i in N+1:2N
            sol .+= coef_vec[i] * sin.(2 * (i - N) * theta)
        end

        return sol
    end

    theta_samples = range(0, 2π, length=500)
    solution = calc_expansion_sol(theta_samples)
    cdf = cumsum(solution)
    jeffery_spline = Spline1D(cdf, theta_samples; k=3)

    # Initialise swimmers.
    function initialize_swimmer!(i, t)
        swim_x[i] = 0.0
        swim_y[i] = sqrt(rand()) * L_y
        swim_theta[i] = jeffery_spline(rand())
        swim_tstart[i] = t
        num_introduced += 1
    end

    num_introduced = 0
    for i in 1:block_size
        initialize_swimmer!(i, 0.0)
    end

    @time begin
    let t_current = 0.0
    while num_introduced < target
        # Simulate block of swimmers
        for i in 1:block_size
        
        
        
            swim_x[i] += (dt * (swim_y[i] + 3*dip*sin(2*swim_theta[i])/(8* swim_y[i]^2)+ Vs * cos(swim_theta[i])))
            swim_y[i] += (dt * (-3*dip.*(1 - 3*(sin(swim_theta[i]))^2)/(8* swim_y[i]^2)+Vs * sin(swim_theta[i])))
            swim_theta[i] += dt * (-3*dip*sin.(2*swim_theta[i]).*(1 .+ beta/2*(1 .+ sin.(swim_theta[i])).^2)./(16* swim_y[i]^3)-0.5 + 0.5 * beta * cos(2 * swim_theta[i])) + D_step * randn()

            if swim_y[i] > L_y
                swim_y[i] = 2 * L_y - swim_y[i]
                swim_theta[i] = -swim_theta[i]
            elseif swim_y[i] < 0
                #if rand() < kappa
                    # Stick to the wall.
                    # Record position and time of sticking
                   push!(stick_events, (swim_x[i], t_current + dt - swim_tstart[i]))
                    # Reinitialize swimmer
                    initialize_swimmer!(i, t_current + dt)
                #else
                    # Reflect off the wall
                   # swim_y[i] = -swim_y[i]
                    #swim_theta[i] = -swim_theta[i]
                end
            end

            if swim_x[i] > L_x
                # Restart swimmer
                initialize_swimmer!(i, t_current + dt)
            end

        end

        t_current += dt

    end
    end
    end
    return stick_events
end

# Do a sweep.
rates = Float64[]

num_params = min(length(gammas), length(Vss), length(Pers))
for i in 1:num_params
    gamma = gammas[i]
    Vs = Vss[i]
    Pe_r = Pers[i]
    dip = dips[i]
    stick_events = run_sim(Pe_r, Vs, gamma, kappa,dip)

    # Calculate the rate of stick events
    rate = gamma * length(stick_events) / (target / cellsPerSecond)
    push!(rates, rate)  

  #  println("--------------------------------")
    # Print the parameters:
   # println("Parameters:")
    #println("Pe_r = $Pe_r")
    #println("Vs = $Vs")
    #println("gamma = $gamma")
    #println("beta = $beta")
    #println("kappa = $kappa")

    # If we assume 10000 swimmers added to the domain per second, this gives a collision rate of
    #p#rintln("Estimated collision rate at $cellsPerSecond cells/s: $rate collisions per second.")
    println(i,"--------------------------------")
end

writedlm("data-rates/Vs$Vs_val-beta-$beta-D_r-$D_r-rep-$rep-rate-dip.txt", rates, "   ")


#println("Summary of results:")
#for i in 1:num_params
    #println("Pe_r=$(Pers[i]), Vs=$(Vss[i]), gamma=$(gammas[i]) => rate=$(rates[i])")
#end