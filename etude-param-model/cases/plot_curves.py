#!/usr/bin/env python3
"""Courbes force-ouverture de tous les cas, disposées en carte dans l'espace (xi, zeta).

À lancer dans le répertoire des cas (celui qui contient cases.txt) :

    python3 plot_curves.py

Une case par cas : xi en colonne (croissant vers la droite), zeta en ligne (croissant vers le haut), comme dans
les cartes de plot_maps.py. Toutes les cases ont les mêmes axes :
  ouverture de 0 à 1.1 x la plus grande ouverture à la rupture (au pic de force) de tous les cas ;
  force     de 0 à 1.1 x la plus grande force de tous les cas.
La force et l'ouverture sont lues comme dans plot_maps.py (top.txt, bottom.txt). Sur chaque courbe, le point
plein marque le pic de force, le point creux le premier lien rompu (breakHistory.txt).

Produit : courbes_force_ouverture.png
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_maps import INK, INK2, GRID, SURFACE, analyse, load_curve, read_cases

SERIES = "#2a78d6"


def main():
    names, cases = read_cases()
    if names is None or len(names) < 2:
        sys.exit("plot_curves.py : cases.txt doit avoir deux paramètres (xi, zeta)")
    px, py = names[0], names[1]

    data = {}
    for d, p in cases:
        curve, r = load_curve(d), analyse(d)
        if curve is None or r is None:
            print(f"  {d} : non terminé ou inexploitable, case vide")
            continue
        data[(p[px], p[py])] = (curve, r)
    if not data:
        sys.exit("plot_curves.py : aucun cas exploitable")

    # mêmes axes pour toutes les cases
    d_max = 1.1 * max(r["D_pic"] for _, r in data.values()) * 1e6           # µm
    f_max = 1.1 * max(float(np.max(c[1])) for c, _ in data.values()) * 1e3  # mN

    xs = sorted({p[px] for _, p in cases})
    ys = sorted({p[py] for _, p in cases}, reverse=True)  # zeta croissant vers le haut
    plt.rcParams.update({"font.size": 9, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})
    fig, axs = plt.subplots(len(ys), len(xs), figsize=(2.3 * len(xs) + 1.2, 1.9 * len(ys) + 1.0),
                            sharex=True, sharey=True, squeeze=False)

    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            ax = axs[j][i]
            ax.set_xlim(0, d_max)
            ax.set_ylim(0, f_max)
            ax.grid(True, color=GRID, linewidth=0.6)
            ax.tick_params(colors=INK2, labelsize=8, length=2)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color("#c9c8c0")
            if (x, y) not in data:
                ax.text(0.5, 0.5, "—", transform=ax.transAxes, ha="center", va="center", color=INK2)
                continue
            (D, F), r = data[(x, y)]
            ax.plot(D * 1e6, F * 1e3, color=SERIES, linewidth=1.6, clip_on=True)
            ax.plot(r["D_pic"] * 1e6, r["F_pic"] * 1e3, "o", color=SERIES, markersize=5,
                    markeredgecolor=SURFACE, markeredgewidth=1.0)
            if np.isfinite(r["D_1re"]) and r["D_1re"] * 1e6 <= d_max:
                ax.plot(r["D_1re"] * 1e6, float(np.interp(r["D_1re"], D, F)) * 1e3, "o", markersize=5,
                        markerfacecolor=SURFACE, markeredgecolor=SERIES, markeredgewidth=1.2)

    for i, x in enumerate(xs):
        axs[-1][i].set_xlabel("ouverture (µm)", color=INK2, fontsize=8)
        axs[0][i].set_title(f"$\\xi$ = {x:g}", color=INK, fontsize=10)
    for j, y in enumerate(ys):
        axs[j][0].set_ylabel(f"$\\zeta$ = {y:g}\n\nforce (mN)", color=INK, fontsize=9)

    handles = [
        plt.Line2D([], [], color=SERIES, linewidth=1.6, label="force – ouverture"),
        plt.Line2D([], [], color=SERIES, marker="o", linestyle="", markersize=5, label="pic de force"),
        plt.Line2D([], [], color=SERIES, marker="o", linestyle="", markersize=5, markerfacecolor=SURFACE,
                   label="premier lien rompu"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.suptitle("Courbes force – ouverture dans l'espace ($\\xi$, $\\zeta$)", color=INK, fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    fig.savefig("courbes_force_ouverture.png", dpi=150)
    print(f"{len(data)} cas -> courbes_force_ouverture.png "
          f"(axes communs : ouverture 0-{d_max:.2f} µm, force 0-{f_max:.2f} mN)")


if __name__ == "__main__":
    main()
