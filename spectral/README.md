# Spectral solution (C++)

Exact three-point return correlations from the Bethe eigenspectrum, extending
[BetheSF](https://data.mendeley.com/datasets/3bs74vf72n/2) by Lapolla and Godec.

**All of my code is in one file — [`BetheSF/SingleFileCorrelations.cpp`](BetheSF/SingleFileCorrelations.cpp).**
Everything else under `BetheSF/` is theirs, except for a fenced block of
declarations in `SingleFile.hpp` (member functions have to be declared inside
the class). See [`../NOTICE.md`](../NOTICE.md) for the full account.

## Build and run

```bash
make && ./example
```

`example` computes `C` on a 40×40 grid of waiting and lag times for `N=5`,
tagged particle 3, return point `x=0.5`, intermediate point `y=0.6`, and writes
`correlation.csv`. Takes about 9 s single-threaded.

```bash
./example [N] [tag] [x] [y] [nbins_t] [max_many_eig] [n_threads]
```

OpenMP is detected automatically. Apple's `clang` (which `g++` aliases to on
macOS) does not support `-fopenmp`, so on a Mac use real GCC to get threading:

```bash
make CXX=g++-14
```

On an 8-thread build the default example goes from 8.7 s to 4.0 s. The
parallelism is only in `precompute_Vkl`.

## What the code does

The spectral representation of the correlation function is

```
C(x,tau,y,t) = Psi00(x)^-2  sum_{k,l}  Psi0k(x) Psil0(x)
                 [ Psikl(y) Psi00(x) - Psi0l(y) Psik0(x) e^{-lambda_k t} ]
                 e^{-lambda_k tau - lambda_l t}
```

where `k`, `l` are multi-indices — one index per particle — and `Psikl` are
overlap elements between many-body eigenfunctions.

Two things make this expensive. The multi-indices mean the truncated sum has
factorially many terms in `N`. And while `Psi0k`-type elements admit
combinatorial shortcuts (which upstream implements), the general `Psikl` do not.

`correlation_looped` addresses this as far as the structure allows:

- `Vkl(x1)` depends on neither `x`, `t` nor `tau`, so the whole matrix is
  precomputed once and reused across the entire grid — this is what
  `precompute_Vkl` does, and it is the part that parallelises
- the remaining overlap elements are hoisted as far out of the loop nest as
  their arguments permit
- `V0k` is obtained from `Vk0` by a multiplicity ratio rather than recomputed

Before any threading, this alone is much faster than calling the scalar
`correlation` in a loop. Measured on a 6×6 grid of times at `max_many_eig=25`,
single-threaded:

```
scalar in a loop  : 5346 ms
correlation_looped:  160 ms      33x, checksums agree to 3e-15
```

The margin grows with grid size, since the `Vkl` precomputation is amortised
over more evaluation points.

### Eigenvalue-shell selection

`correlation_looped` takes an `eig_to_exclude` argument, which drops all `(k,l)`
pairs where **both** indices are in the given set. Passing the eigenfunctions of
a smaller `max_many_eig` therefore returns the *incremental* contribution of the
higher shell, so that

```
full sum  =  run(max_many_eig = K, exclude = eigenfunctions of k)
           + run(max_many_eig = k, exclude = {})
```

This lets a convergence study be split across independent jobs and summed.

## Limitations

Being honest about where this stands — this is thesis code from 2021, tidied but
not rewritten.

- **It does not scale.** `Psikl` cannot be reduced below factorial complexity in
  `N`, and in practice `C` is not computable much beyond `N = 5`. Measured
  wall-clock on the group cluster for one convergence run (30 return points,
  30 x 30 times, one intermediate point):

  | N | max_many_eig | runtime |
  |---|---|---|
  | 5 | 180 | 94 h |
  | 6 | 200 | 168 h |
  | 7 | 160 | 446 h |
  | 8 | 110 | 349 h |
  | 9 | 90 | 478 h |

  This is *why* the Brownian dynamics code in [`../bd`](../bd) exists. None of
  the figures in the thesis came from this program.

- **Truncation error at short times.** The spectral sum must be cut off at
  `max_many_eig`; below roughly `t ~ 1/(D pi^2 kmax^2)` the result develops
  small oscillations and can go slightly negative where it should be zero.
  Reaching `t = 5e-5` needed `max_many_eig = 400`. Always check convergence by
  increasing it.

- **Threading is limited to `precompute_Vkl`**, and uses static scheduling even
  though per-iteration cost varies a lot with the multiset. Measured on the
  cluster, going from 50 to 100 threads made things *slower*
  (61,030 s → 82,480 s at `max_many_eig = 400`), which is consistent with
  oversubscription plus load imbalance. `schedule(dynamic)` would likely help;
  I never tried it.

- **Compiler warnings.** The vendored code emits a handful of
  `-Wdefaulted-function-deleted` warnings, because `SingleFile` has a `const`
  member and so its assignment operators cannot be defaulted. These are upstream
  and harmless; I have left them rather than patch someone else's code.

- Only `SingleFileFlat` (free particles in a box) is exercised here.
  `SingleFileOnSlope` is inherited from upstream and untested by me;
  `SingleFileHarmonic` needs `std::hermite` and is unavailable on libc++.
