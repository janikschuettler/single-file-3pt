# Cluster jobs

These are the **actual job scripts as submitted** to the Sun Grid Engine cluster
of the Max Planck Institute for Biophysical Chemistry (now MPI for
Multidisciplinary Sciences) in 2021, kept verbatim. They are included
as a record of how the results were produced, not as a supported entry point —
the drivers they invoke (`main.jl`, `g2p.jl`, `nm`, `conv`) were the untrimmed
versions of what is now `bd/example.jl` and `spectral/example.cpp`. Do not
expect them to run as-is outside that cluster.

| Script | What it ran | Array | Walltime |
|---|---|---|---|
| `bd_correlations.sh` | BD three-point correlations, sweeping `(N, tag, x, deviation)` | `-t 1-25` | 80 h |
| `bd_propagators.sh` | BD two-point propagators, Markovian and not | `-t 1-6` | 50 h |
| `spectral_correlations.sh` | spectral `C`, sweeping the intermediate point `x1` | `-t 1-21` | 80 h |
| `spectral_convergence.sh` | spectral convergence in `max_many_eig` | `-t 1-11` | **500 h** |

## How the code was shaped by this

The two least obvious design decisions in this repository both exist because the
work ran as cluster arrays rather than as single jobs.

**Integer accumulators (`bd/`).** `SFBD.bd` returns unnormalised integer counts
and defers all normalisation to `SFBD.correlation`. Histograms of integers add
*exactly*, so the 25 tasks of an array job can each write their counts into a
shared HDF5 file, and the result is bit-for-bit what one enormous run would have
given — no floating-point drift over tens of millions of trajectories. This is
what `bd/accumulate.jl` implements, and it is why a target of 16,000,000
trajectories is reachable at all: it is not one job, it is many, resumable and
extendable at any time by submitting more.

**Eigenvalue-shell decomposition (`spectral/`).** `correlation_looped` takes an
`eig_to_exclude` set and skips `(k,l)` pairs where *both* indices fall in it.
Excluding a low shell therefore returns the incremental contribution of the
higher one, so

```
full sum  =  run(max_many_eig = K, exclude = shell k)  +  run(max_many_eig = k)
```

A convergence study that would otherwise be one un-checkpointable multi-week job
splits into independent tasks that sum afterwards — which matters when the queue
limit is 500 hours and a single `N=9` point takes 478 of them (nearly three weeks).

## Scale

Trajectory counts behind the published figures (thesis Table 1):

| Model / initial condition | N=5 | N=11 | N=21 | N=31 | N=61 |
|---|---|---|---|---|---|
| Markovian | 16,000,000 | ~11,000,000 | 6,000,000 | 1,000,000 | 2,000,000 |
| non-Markovian, equilibrium | 7,000,000 | 1,000,000 | 3,000,000 | 500,000 | 500,000 |
| non-Markovian, non-equilibrium | 2,000,000 | 1,000,000 | 3,000,000 | 1,500,000 | 500,000 |

About 56 million trajectories in total. Each is 100,000 Euler–Maruyama steps
(`T = 5`, `dt = 5e-5`) for every particle in the file, which comes to roughly
**3.5 x 10^13 particle-position updates**.

Wall-clock for the spectral solver, which is the reason the simulations exist at
all (one convergence run: 30 return points, 30 x 30 times):

| N | max_many_eig | runtime |
|---|---|---|
| 5 | 180 | 94 h |
| 6 | 200 | 168 h |
| 7 | 160 | 446 h |
| 8 | 110 | 349 h |
| 9 | 90 | 478 h |

Thread scaling of `precompute_Vkl` on that hardware was poor beyond a point:
at `max_many_eig = 400`, 50 threads took 61,030 s and 100 threads took
82,480 s. Static scheduling over iterations whose cost varies with the multiset
is the likely cause; see [`../spectral/README.md`](../spectral/README.md).
