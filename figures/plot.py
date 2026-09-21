#!/usr/bin/env python3
"""
Builds figure.png from the CSVs produced by the two examples.

    spectral/example  ->  spectral_nm.csv
    bd/example.jl nm  ->  bd_nm.csv
    bd/example.jl m   ->  bd_m.csv

Left panel:  C(t,tau) from the exact spectral solution.
Right panel: a cut at fixed t, comparing the spectral solution with Brownian
             dynamics for the same (non-Markovian) system, and with Brownian
             dynamics for the Markovian approximation.

Usage: python3 plot.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
CUT_T = 0.1          # waiting time for the right-hand panel (~ the relaxation time)


def load(name):
    """Read t, tau, C from one CSV into a grid."""
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("t,"):
                continue
            f = line.split(",")
            rows.append((float(f[0]), float(f[1]), float(f[2])))
    a = np.array(rows)
    t, tau = np.unique(a[:, 0]), np.unique(a[:, 1])
    C = np.full((t.size, tau.size), np.nan)
    ti = {v: i for i, v in enumerate(t)}
    tj = {v: j for j, v in enumerate(tau)}
    for r in a:
        C[ti[r[0]], tj[r[1]]] = r[2]
    return t, tau, C


def row_at(t, tau, C, target):
    """The row of C nearest to t = target."""
    i = int(np.argmin(np.abs(np.log(t) - np.log(target))))
    return t[i], tau, C[i]


spectral = load("spectral_nm.csv")
bd_nm = load("bd_nm.csv")
bd_m = load("bd_m.csv")

if spectral is None:
    raise SystemExit("spectral_nm.csv not found -- run spectral/example first")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

# ---- left: the exact correlation map -------------------------------------
t, tau, C = spectral
lim = np.nanpercentile(np.abs(C), 99)
mesh = ax1.pcolormesh(tau, t, C, cmap="RdBu_r", shading="nearest",
                      norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim))
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel(r"lag time  $\tau$")
ax1.set_ylabel(r"waiting time  $t$")
ax1.set_title("exact (spectral solution)")
fig.colorbar(mesh, ax=ax1, label=r"$C(x,\tau,y,t)$")

t_cut = row_at(*spectral, CUT_T)[0]
ax1.axhline(t_cut, color="k", lw=0.8, ls="--")

# ---- right: a cut, both methods and both models --------------------------
_, tau_s, C_s = row_at(*spectral, CUT_T)
ax2.plot(tau_s, C_s, "k-", lw=1.5, label="non-Markovian, spectral (exact)")

if bd_nm is not None:
    _, tau_b, C_b = row_at(*bd_nm, CUT_T)
    ax2.plot(tau_b, C_b, "o", ms=4, color="C0", alpha=0.8,
             label="non-Markovian, Brownian dynamics")
if bd_m is not None:
    _, tau_b, C_b = row_at(*bd_m, CUT_T)
    ax2.plot(tau_b, C_b, "s", ms=4, color="C3", alpha=0.8,
             label="Markovian approximation, Brownian dynamics")

ax2.axhline(0, color="0.6", lw=0.8)
ax2.set_xscale("log")
ax2.set_xlabel(r"lag time  $\tau$")
ax2.set_ylabel(r"$C(x,\tau,y,t)$")
ax2.set_title(rf"cut at $t = {t_cut:.3g}$")
ax2.legend(frameon=False, fontsize=8)

fig.suptitle(r"Three-point return correlations, $N=5$, tag 3, $x=0.5$, $y=0.6$",
             fontsize=11)
fig.tight_layout()
out = os.path.join(HERE, "figure.png")
fig.savefig(out, dpi=150)
print("wrote", out)
