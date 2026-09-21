"""
* 2021 Janik Schüttler
* Brownian Dynamics simulations for the single-file model. Computes multi-point
* correlation functions.
"""

using LinearAlgebra
using Random
using Statistics

Random.seed!(5)

module SFBD

include("sfbd_utilities.jl")

"""
Defines the Single File type.
"""
struct SingleFile
 N::Int
 tag::Int
 x0::Float64

 ic_confinement::Float64
 force_const::Float64
end

# struct SingleFile
#  N::Int
#  tag::Int
#  x0::Float64
#  parameter::Dict
# end
# SingleFile(N::Int, tag::Int, x0::Float64) = SingleFile(N,tag,x0,Dict())

"""
Defines the BDSimulation type that contains all necessary information to perform
a general Brownian dynamics simulation.
"""
struct BDSimulation{T1, T2, T3}
 N::Int
 force!::T1
 diffusion!::T2
 sample_initial_condition!::T3
end

"""
Force and diffusion functions used in the simulations.
- force_zero!: a force that is constant zero
- force_markov!: the Markovian force that ensures the same steady-state
                 distribution as for the tagged-particle dynamics
- diffusion_unit!: a diffusion that is constant one
"""
force_zero!(force,x,t,sf) = for i in 1:length(force) force[i] += 0.0 end
force_const!(force,x,t,sf) = for i in 1:length(force) force[i] += -sf.force_const end
force_markov!(force,x,t,sf) = for i in 1:length(force) force[i] += (sf.tag-1)/(x[i]+1e-10) - (sf.N-sf.tag)/(1-x[i]-1e-10) end
function force_const_and_markov!(force,x,t,sf)
  force_const!(force,x,t,sf)
  force_markov!(force,x,t,sf)
end
diffusion_unit!(force,x,t,sf) = for i in 1:length(force) force[i] += 1.0 end

"""
Samples equilibrium initial conditions for all N particle in the single-file such
that the tagged particle is pinned at sf.x0 and all other particles are sampled
uniformly given that the tagged particle is fixed, i.e. such that tag-1
particles are to the left of tag and N-tag to the right of the tagged particle.
"""
function projected_initial_condition!(pos, tag, x0)
  for particle in 1:length(pos)
    if particle < tag
      pos[particle] = x0 * rand()
    elseif particle == tag
      pos[particle] = x0
    else
      pos[particle] = x0 + (1-x0) * rand()
    end
  end
  sort!(pos)
end
projected_initial_condition!(pos, sf) = projected_initial_condition!(pos, sf.tag, sf.x0)

"""
Samples non-equilibrium initial conditions for all N particles in the single-file
such that the tagged particle is pinned at sf.x0 and all other particles are
sampled uniformly in a neighborhood around sf.x0.
"""
function confined_initial_condition!(pos, tag, x0, dx)
  boundary_l, boundary_r = x0 - dx, x0 + dx
  boundary_l = boundary_l < 0.0 ? 0.0 : boundary_l
  boundary_r = boundary_r > 1.0 ? 1.0 : boundary_r
  for particle in 1:length(pos)
    if particle < tag
      pos[particle] = boundary_l + (x0 - boundary_l) * rand()
    elseif particle == tag
      pos[particle] = x0
    else
      pos[particle] = x0 + (boundary_r - x0) * rand()
    end
  end
  sort!(pos)
end
confined_initial_condition!(pos, sf) = confined_initial_condition!(pos, sf.tag, sf.x0, sf.ic_confinement)

"""
Returns for a tagged-particle of index tag within a single-file of N particle
the steady state probability density at point q.
"""
function steady_state(q, N::Int, tag::Int)
  density = q.^(tag-1) .* (1 .- q).^(N-tag)
  if N <= 20
    density *= factorial(N) / factorial(tag-1) / factorial(N-tag)
  else
    larger, smaller = max(tag-1, N-tag), min(tag-1, N-tag)
    for i in 1:smaller
      density *= (larger+1+i) / i
    end
    density *= larger+1
  end
  return density
end
steady_state(q; sf) = steady_state(q, sf.N, sf.tag)

"""
Performs a Brownian dynamics simulation of a Single-File model.
"""
function bd(; sf, bdsimulation, M, time, space, time_coarsening_indices, x1)
  trajectory, trajectory_bin = zeros(Float64, bdsimulation.N, length(time)), zeros(Int, length(time))
  tmp_position = zeros(Float64, bdsimulation.N)

  force, diffusion = zeros(Float64, bdsimulation.N), ones(Float64, bdsimulation.N)

  x_bin, x1_bins = bin(sf.x0, space), zeros(Int, length(x1))
  increment_x = x_bin == 1 ||  x_bin == length(space) ? Int(2/step(space)) : Int(1/step(space))
  increment_x1s = zeros(Float64, length(x1))
  for i in 1:length(x1)
    x1_bins[i] = bin(x1[i], space)
    increment_x1s[i]  = x1_bins[i] == 1 || x1_bins[i] == length(space) ? Int(2/step(space)) : Int(1/step(space))
  end

  time_coarsening_indices_sums = zeros(Int, length(time_coarsening_indices)^2+length(time_coarsening_indices)) .- 1
  all_sums!(time_coarsening_indices_sums, time_coarsening_indices)

  density2p = zeros(Int, 1+length(x1), length(time_coarsening_indices_sums))
  g3ps = zeros(Int, length(x1), length(time_coarsening_indices), length(time_coarsening_indices))
  # correlations, g2ps = zeros(size(g3ps)), zeros(size(g3ps))

  # Threads.@threads
  for _ in 1:M
    for (i_t, t) in enumerate(time)
        # handle initial condition
        if i_t == 1
          bdsimulation.sample_initial_condition!(tmp_position, sf)
          for particle in 1:bdsimulation.N
            trajectory[particle, 1] = tmp_position[particle]
            if particle == sf.tag || bdsimulation.N == 1
              trajectory_bin[i_t] = bin(trajectory[particle, i_t], space)
            end
          end
          continue
        end

        # update force & diffusion vectors
        fill!(force, 0.0)
        bdsimulation.force!(force, tmp_position, t, sf)
        fill!(diffusion, 0.0)
        bdsimulation.diffusion!(diffusion, tmp_position, t, sf)

        for particle in 1:bdsimulation.N
          # integration step
          tmp_position[particle]  = trajectory[particle, i_t-1]
          tmp_position[particle] += step(time) * diffusion[particle] * force[particle]
          tmp_position[particle] += sqrt(2*step(time)*diffusion[particle]) * randn()

          # boundary condition: no flux at x=0.0 and x=1.0
          if tmp_position[particle] < 0.0 || tmp_position[particle] > 1.0
            tmp_position[particle] = trajectory[particle, i_t-1]
          end
        end

        # no flux boundary condition within particles
        sort!(tmp_position)

        for particle in 1:bdsimulation.N
          trajectory[particle, i_t] = tmp_position[particle]
          if particle == sf.tag || bdsimulation.N == 1
            trajectory_bin[i_t] = bin(trajectory[particle, i_t], space)
          end
        end
    end

    for (i, i_t) in enumerate(time_coarsening_indices_sums)
      if i_t > 0
        i_pos = trajectory_bin[i_t]
        if i_pos == x_bin
          density2p[1, i] += increment_x
        end
        for (i_x1, x1_bin) in enumerate(x1_bins)
          if i_pos == x1_bin
            density2p[1+i_x1, i] += increment_x1s[i_x1]
          end
        end
      end
    end

    for (i_t1_coarse, i_t1) in enumerate(time_coarsening_indices)
      if i_t1 == 1 continue end
      window_size = 0 #Int(round(5 * i_t / length(tiime)))
      for (i_x1, x1_bin) in enumerate(x1_bins)
        # if trajectory_bin[i_t1] == x1_bin
        if is_in(trajectory_bin, i_t1, x1_bin, window_size)
          for (i_tau_coarse, i_tau) in enumerate(time_coarsening_indices)
            if i_t1+i_tau > length(time)
              break
            end
            # if trajectory_bin[i_t1+i_tau] == x_bin
            if is_in(trajectory_bin, i_t1+i_tau, x_bin, window_size)
              g3ps[i_x1, i_t1_coarse, i_tau_coarse] += increment_x * increment_x1s[i_x1]
            end
          end
        end
      end
    end
  end

  g3ps, density2p, trajectory
end


function bd_g2p(; sf, bdsimulation, M, time, space, time_coarsening_indices)
  trajectory, trajectory_bin = zeros(Float64, bdsimulation.N, length(time)), zeros(Int, length(time))
  tmp_position = zeros(Float64, bdsimulation.N)

  force, diffusion = zeros(Float64, bdsimulation.N), ones(Float64, bdsimulation.N)

  increment = Int(1 / step(space))
  density2p = zeros(Int, length(space), length(time_coarsening_indices), length(space))

  # Threads.@threads
  for (i_x0, x0) in enumerate(space)
    for _ in 1:M
      for (i_t, t) in enumerate(time)
          # handle initial condition
          if i_t == 1
            if bdsimulation.N == 1
              confined_initial_condition!(tmp_position, 1, x0, sf.ic_confinement)
            else
              confined_initial_condition!(tmp_position, sf.tag, x0, sf.ic_confinement)
            end
            for particle in 1:bdsimulation.N
              trajectory[particle, 1] = tmp_position[particle]
              if particle == sf.tag || bdsimulation.N == 1
                trajectory_bin[i_t] = bin(trajectory[particle, i_t], space)
              end
            end
            continue
          end

          # update force & diffusion vectors
          fill!(force, 0.0)
          bdsimulation.force!(force, tmp_position, t, sf)
          fill!(diffusion, 0.0)
          bdsimulation.diffusion!(diffusion, tmp_position, t, sf)

          for particle in 1:bdsimulation.N
            # integration step
            tmp_position[particle]  = trajectory[particle, i_t-1]
            tmp_position[particle] += step(time) * diffusion[particle] * force[particle]
            tmp_position[particle] += sqrt(2*step(time)*diffusion[particle]) * randn()

            # boundary condition: no flux at x=0.0 and x=1.0
            if tmp_position[particle] < 0.0 || tmp_position[particle] > 1.0
              tmp_position[particle] = trajectory[particle, i_t-1]
            end
          end

          # no flux boundary condition within particles
          sort!(tmp_position)

          for particle in 1:bdsimulation.N
            trajectory[particle, i_t] = tmp_position[particle]
            if particle == sf.tag || bdsimulation.N == 1
              trajectory_bin[i_t] = bin(trajectory[particle, i_t], space)
            end
          end
      end

      for (i_t, t) in enumerate(time_coarsening_indices)
        i_x = trajectory_bin[time_coarsening_indices[i_t]]
        if i_x == 1 || i_x == length(space)
          density2p[i_x0, i_t, i_x] += 2*increment
        else
          density2p[i_x0, i_t, i_x] += increment
        end
      end
    end
  end

  density2p, trajectory
end


end
