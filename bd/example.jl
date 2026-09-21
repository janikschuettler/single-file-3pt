#=
Minimal example: three-point return correlations from Brownian dynamics.

Simulates M independent trajectories of a single file of N hard-core particles
on [0,1], tracks the tagged particle, and estimates

    C(x,tau,y,t) = Q(x,t+tau; y,t | x) - Q(x,t+tau|x) Q(y,t|x)

Two models share one integrator, selected by `model`:
  :nm  the full N-particle single file  -- tagged-particle dynamics is
       non-Markovian, because the other N-1 particles have been projected out
  :m   one particle in the potential of mean force U_M = -log(P_ss), a
       Markovian process with the *same* stationary distribution

Writes correlation_bd_<model>.csv.

Usage:  julia example.jl [model] [N] [tag] [x] [y] [M]
Default: julia example.jl nm 5 3 0.5 0.6 200000
=#

include("sfbd.jl")

using Random
using Printf

model = length(ARGS) > 0 ? Symbol(ARGS[1]) : :nm
N     = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 5
tag   = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 3
x     = length(ARGS) > 3 ? parse(Float64, ARGS[4]) : 0.5
y     = length(ARGS) > 4 ? parse(Float64, ARGS[5]) : 0.6
M     = length(ARGS) > 5 ? parse(Int, ARGS[6]) : 200_000

Random.seed!(5)

# Discretisation. dt and the space bin width are the two systematic errors here;
# everything else is Monte Carlo noise that shrinks as 1/sqrt(M).
space   = range(0, 1, length = 51)
time_bd = range(0, 2, step = 1e-4)
t_bd_to_plt = SFBD.coarsen_time(time_bd, 40, 10, fraction = 0.0)
time_plt    = time_bd[t_bd_to_plt]

# confinement = 1.0 means the equilibrium initial condition: the tagged particle
# is pinned at x and the remaining particles are drawn uniformly on either side
# of it. Smaller values start the file in a non-equilibrium, confined state.
sf = SFBD.SingleFile(N, tag, x, 1.0, 0.0)

simulation = if model == :nm
    SFBD.BDSimulation(N, SFBD.force_const!, SFBD.diffusion_unit!,
                      SFBD.confined_initial_condition!)
elseif model == :m
    SFBD.BDSimulation(1, SFBD.force_const_and_markov!, SFBD.diffusion_unit!,
                      (pos, sf) -> pos[1] = sf.x0)
else
    error("model must be :nm or :m, got :$model")
end

@printf("model = %s, N = %d, tag = %d, x = %.3f, y = %.3f, M = %d\n",
        model, N, tag, x, y, M)
@printf("%d integration steps of dt = %g, %d sampled times\n",
        length(time_bd), step(time_bd), length(t_bd_to_plt))

kwargs = (bdsimulation = simulation, sf = sf, time = time_bd, space = space,
          x1 = [y], time_coarsening_indices = t_bd_to_plt)

SFBD.bd(; kwargs..., M = 2)                      # warm up the JIT
elapsed = @elapsed out = SFBD.bd(; kwargs..., M = M)
@printf("simulated %d trajectories in %.1f s\n", M, elapsed)

# bd returns *unnormalised* integer counts, so results from independent runs can
# be summed exactly before normalising. correlation() does the normalisation and
# assembles the connected and disconnected parts.
connected_nn, propagators_nn = out[1], out[2]
C, connected, disconnected =
    SFBD.correlation(connected_nn[1, :, :], propagators_nn[[1, 2], :],
                     M, t_bd_to_plt)

open("correlation_bd_$(model).csv", "w") do io
    println(io, "# model=$model N=$N tag=$tag x=$x y=$y M=$M dt=$(step(time_bd))")
    println(io, "t,tau,C,connected,disconnected")
    for i in eachindex(time_plt), j in eachindex(time_plt)
        @printf(io, "%g,%g,%g,%g,%g\n", time_plt[i], time_plt[j],
                C[i, j], connected[i, j], disconnected[i, j])
    end
end

println("wrote correlation_bd_$(model).csv")
