#!/usr/bin/env python3
"""Figures propres au compte rendu (à lancer dans etude-param-frank/rapport après les scripts de cases/).

  fig_raideur_serie.png   raideur mesurée K / Kcell en fonction de zeta, comparée à 1/K = 1/Kcell + b/kc
  fig_vitesse_onde.png    vitesse de la fissure principale rapportée à c0 = sqrt(E/rho), en fonction de zeta, et part
                          des ramifications en fonction de v/c0

Lit ../cases/resultats.csv et ../cases/crackwatch.csv ; écrit aussi valeurs_rapport.txt (constantes citées dans le texte).
"""

import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../cases")
from plot_maps import GRID, INK, INK2, RAMP, SURFACE  # noqa: E402

from matplotlib.colors import LogNorm  # noqa: E402

KS = 2000.0                                  # N/m
W, H = 13.6e-3, 16e-3                        # dimensions de l'échantillon (m)
RHO = 13084 * 8e-8 / (W * H)                 # masse surfacique (kg/m²) : 13084 noeuds de 8e-8 kg

rows = {}
for f in ("../cases/resultats.csv", "../cases/crackwatch.csv"):
    for r in csv.DictReader(open(f)):
        rows.setdefault((float(r["XI"]), float(r["ZETA"])), {}).update(r)
keys = sorted(rows)
X = np.array([k[0] for k in keys])
Z = np.array([k[1] for k in keys])
K = np.array([float(rows[k]["K_N_par_m"]) for k in keys])
V = np.array([float(rows[k]["vitesse_m_par_s"]) for k in keys])
R = np.array([float(rows[k]["part_ramifications"]) for k in keys])
Kcell = 12 * X * KS / (1 + 9 * X)

# 1/K = a/Kcell + b/kc, kc = zeta Kcell
(a, b), *_ = np.linalg.lstsq(np.column_stack([1 / Kcell, 1 / (Z * Kcell)]), 1 / K, rcond=None)
err = np.abs(1 / (a / Kcell + b / (Z * Kcell)) / K - 1)
c0 = np.sqrt(K * H / W / RHO)
vr = V / c0

cmap = RAMP
norm = LogNorm(0.1, 10)
plt.rcParams.update({"font.size": 10, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2)


fig, ax = plt.subplots(figsize=(6.2, 4.2))
zz = np.logspace(-1.2, 1.1, 200)
ax.plot(zz, zz / (a * zz + b), color=INK2, linewidth=1.2, linestyle="--",
        label=f"$K/K_{{cell}} = \\zeta / ({a:.3f}\\,\\zeta + {b:.3f})$")
sc = ax.scatter(Z, K / Kcell, c=X, cmap=cmap, norm=norm, s=22, zorder=3, edgecolor="none")
ax.set_xscale("log")
ax.set_xlabel("$\\zeta$")
ax.set_ylabel("$K / K_{cell}$")
ax.legend(frameon=False, loc="lower right")
cb = fig.colorbar(sc, ax=ax)
cb.set_label("$\\xi$")
cb.outline.set_visible(False)
style(ax)
fig.tight_layout()
fig.savefig("fig_raideur_serie.png", dpi=200)
plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
sc = axs[0].scatter(Z, vr, c=X, cmap=cmap, norm=norm, s=22, edgecolor="none", zorder=3)
axs[0].set_xscale("log")
axs[0].set_xlabel("$\\zeta$")
axs[0].set_ylabel("$v / c_0$")
axs[0].set_title("Vitesse de la fissure principale", color=INK)
axs[1].scatter(vr, R, c=X, cmap=cmap, norm=norm, s=22, edgecolor="none", zorder=3)
axs[1].set_xlabel("$v / c_0$")
axs[1].set_ylabel("part des ramifications")
axs[1].set_title("Ramifications et vitesse", color=INK)
for ax in axs:
    style(ax)
cb = fig.colorbar(sc, ax=axs)
cb.set_label("$\\xi$")
cb.outline.set_visible(False)
fig.savefig("fig_vitesse_onde.png", dpi=200, bbox_inches="tight")
plt.close(fig)

with open("valeurs_rapport.txt", "w") as f:
    f.write(f"a = {a:.4f}\nb = {b:.4f}\nerreur max = {100 * err.max():.2f} %\nerreur mediane = {100 * np.median(err):.2f} %\n")
    f.write(f"rho = {RHO:.3f} kg/m2\nc0 min = {c0.min():.2f} m/s, max = {c0.max():.2f} m/s\n")
    for lo, hi in ((0, 1.01), (1.99, 11)):
        m = (Z > lo) & (Z < hi)
        f.write(f"zeta in ]{lo},{hi}[ : v/c0 {vr[m].min():.2f}-{vr[m].max():.2f}, ramif mediane {np.median(R[m]):.2f}\n")
print(open("valeurs_rapport.txt").read())
