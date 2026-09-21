# Attribution

## BetheSF

`spectral/BetheSF/` is **not my code**. It is the reference implementation
published with

> A. Lapolla and A. Godec, *BetheSF: Efficient computation of the exact
> tagged-particle propagator in single-file systems via the Bethe eigenspectrum*
> — Mendeley Data, version 2.
> <https://data.mendeley.com/datasets/3bs74vf72n/2>

Copyright (C) 2020 Alessio Lapolla and Aljaž Godec, MIT licensed. The original
`LICENSE.txt` and `README.txt` are kept in that directory unmodified.

BetheSF diagonalises the Fokker–Planck operator of the full N-particle single
file via a coordinate Bethe ansatz and projects the result onto the coordinate
of one tagged particle. It implements **two-point** propagators.

## What I added

I extended it to compute **three-point return correlations**.

**Every line I wrote lives in one file:
[`spectral/BetheSF/SingleFileCorrelations.cpp`](spectral/BetheSF/SingleFileCorrelations.cpp).**
It carries my copyright header and contains nothing of theirs.

| Function | What it does |
|---|---|
| `joint3dens` | three-point joint density |
| `correlation` (scalar) | `C` at one `(x, t, y, tau)` — the direct transcription of the spectral sum |
| `correlation` (2 vectorised overloads) | the same over grids of times / positions, hoisting the overlap elements out of the grid loops |
| `correlation_looped` | the production version: return correlations over a grid, with precomputed `Vkl` and eigenvalue-shell selection |
| `precompute_Vkl` | parallel (OpenMP) precomputation of the expensive `Vkl(x1)` overlap elements |
| `joint2dens` (vector overload) | propagator over a vector of times |

### The one unavoidable exception

These are **member functions of `SingleFile`**, so C++ requires their
*declarations* to appear inside the class definition — which lives in their
`SingleFile.hpp`. Those declarations are fenced off:

```cpp
//=============================================================================
// Extension by Janik Schuettler (2021), MIT licensed -- three-point return
// correlations. ... Nothing below this marker, down to the matching end
// marker, is part of upstream BetheSF.
//=============================================================================
   ...
//===================== end of Janik Schuettler's extension ===================
```

plus two small accessors (`eigenfunctions()`, `print_eigenfunctions()`) marked
the same way further down. C++ does permit a class's members to be *defined*
across several translation units, which is what makes the separation above
possible at all.

Apart from those marked blocks, `SingleFile.hpp` and every other file under
`spectral/BetheSF/` is upstream code (modulo re-indentation).

### Verifying the boundary

To see exactly which lines are mine, download version 2 of the dataset above and
diff against it:

```bash
diff -uw -B path/to/upstream/SingleFile.cpp spectral/BetheSF/SingleFile.cpp
diff -uw -B path/to/upstream/SingleFile.hpp spectral/BetheSF/SingleFile.hpp
```

`-w -B` is worth using because upstream indents with 4 spaces and this copy uses
2. `SingleFileCorrelations.cpp` has no upstream counterpart — it is entirely new.

Everything under `bd/` is mine.

## Licence

This repository is MIT licensed, the same as the code it builds on.
