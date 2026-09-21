# Brownian dynamics (Julia)

Monte Carlo estimate of three-point return correlations, for the full
single file and for its Markovian approximation.

This is the code that produced the results in the thesis. It is slower per data
point than the [spectral solution](../spectral) but it scales: the thesis used
it up to `N = 61` and 16 million trajectories, where the exact route stops
around `N = 5`.

## Run

The engine (`sfbd.jl`, `sfbd_utilities.jl`) has no dependencies beyond the Julia
standard library. Only `accumulate.jl` needs a package (HDF5.jl); `example.jl`
writes plain CSV and does not.

```bash
julia example.jl                       # non-Markovian, N=5, tag=3, M=200000
julia example.jl m 5 3 0.5 0.6 200000  # the Markovian approximation
```

```bash
julia example.jl [model] [N] [tag] [x] [y] [M]
```

Writes `correlation_bd_<model>.csv`. Expect a few minutes at the default `M`;
see *Limitations* on how many trajectories you actually need.

## The two models

Both run through the same integrator, `SFBD.bd`, which takes force, diffusion
and initial-condition sampler as fields of a `BDSimulation`. Swapping those
three closures is the only difference between:

- **`:nm`** — the full `N`-particle single file. Particles are integrated
  independently, then `sort!`ed at every step, which is exactly equivalent to
  hard-core reflection for indistinguishable walkers. Watching one particle of
  this means projecting away the other `N-1`, and the resulting tagged-particle
  process is non-Markovian.

- **`:m`** — a *single* particle in the potential of mean force
  `U_M = -log P_ss`, giving the force `F_M = (k-1)/q - (N-k)/(1-q)`. This is a
  genuine Markov process constructed to have **exactly the same stationary
  distribution** as the tagged particle above. It is the controlled comparison:
  same equilibrium, no memory.

Keeping both behind one integrator is deliberate — it means any difference
between them is the physics, not the numerics.

## Notes on the implementation

A few choices that are not obvious from reading the code:

- **`bd` returns unnormalised integer counts**, and normalisation is deferred to
  `SFBD.correlation`. This is what makes results from independent runs summable
  *exactly* — no floating-point drift over millions of trajectories, and
  separate cluster jobs can be added together and normalised once at the end.
  `accumulate.jl` is that merge layer: it adds a run's counts into whatever an
  HDF5 file already holds, refusing to merge runs whose discretisation differs.
  It is how the 16,000,000-trajectory results in the thesis were reached — not
  as one job, but as resumable array tasks. See [`../cluster`](../cluster).

- **`all_sums!` is the reason the disconnected part is unbiased.** The
  disconnected term needs the propagator at `t + tau`, but `t` and `tau` both
  live on a log-spaced coarse grid, so `t + tau` generally is *not* a grid
  point. Rather than interpolating, the code precomputes the index set
  `{t_i} ∪ {t_i + t_j}` and accumulates the single-time propagator at exactly
  those extra indices. Both terms of `C` then come off the same trajectories at
  the exact right times.

- **Julia rather than C++ was a measured choice, not a default.** I wrote an
  equivalent C++ integrator and benchmarked the two; the Julia version was at
  least as fast, given care over allocation. The inner loop allocates nothing:
  buffers are preallocated outside the trajectory loop and mutated in place,
  and observables are accumulated only at the sampled times rather than every
  step.

- **Only the tagged particle's bin index is retained** for the observables; the
  full trajectory matrix is kept only so the previous step is available.

- Time is coarsened logarithmically (`coarsen_time`) because the physics spans
  decades: collisions at `~1/(D N^2)`, relaxation at `~1/(D pi^2)`.

## Limitations

Honest notes — this is thesis code from 2021, tidied but not rewritten.

- **`C` is a difference of two estimated quantities**, so statistical error
  propagates and amplifies. Getting smooth correlation maps takes a lot of
  trajectories. The thesis used, for `N = 5`: 16,000,000 (Markovian) and
  7,000,000 (non-Markovian, equilibrium initial condition). The default `M` in
  `example.jl` is far below that and will look noisy — it is sized to finish in
  minutes, not to be publication quality.

- **`bin()` is a linear scan** over the space grid, even though that grid is
  uniform and an O(1) index would do. It sits in the hot loop. Replacing it with
  `round(Int, (x - first(grid))/step(grid)) + 1` measures ~28% faster overall on
  identical output. I have deliberately left it as it was.

- **Not parallelised.** There is a commented-out `Threads.@threads` on the
  trajectory loop, which is embarrassingly parallel. The thesis parallelised
  across cluster array jobs instead and relied on the exact summability above,
  so this was never needed — but it is the obvious first speedup.

- **Fixed Euler–Maruyama with step rejection at the walls.** A step that would
  leave `[0,1]` is rejected rather than reflected, which is a slight bias at
  finite `dt`. The Markovian force is regularised with `eps = 1e-10` since it
  diverges at the boundaries for interior particles.

- `bd_g2p` loops over all initial positions in the space grid, which is
  wasteful if you only want one. It is kept because it is what generated the
  propagator data.

- No test suite. The one correctness check that exists is the comparison against
  the spectral solution described in the [top-level README](../README.md), which
  I ran while assembling this repository — not during the thesis.
