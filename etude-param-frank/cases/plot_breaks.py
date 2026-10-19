#!/usr/bin/env python3
"""Cartes et courbes de la rupture (breakEvol.txt, breakHistory.txt) dans l'espace (xi, zeta).

À lancer dans le répertoire des cas (celui qui contient cases.txt) :

    python3 plot_breaks.py [--ligament 12.42e-3]

Pour chaque cas terminé :
  breakHistory.txt  une ligne par face collée rompue : instant, cellules et noeuds, released_NRJ = W, énergie élastique
                    stockée dans les deux points de colle de la face au moment de la rupture (J)
  breakEvol.txt     toutes les nstepPeriodRecord pas : instant, somme des taux de restitution G des faces rompues
                    (J/m, non utilisée ici), longueur cumulée des faces rompues (m)

Grandeurs :
  longueur rompue / ligament   longueur finale des faces rompues divisée par la longueur d'une fissure droite de la
                               pointe de l'entaille au bord opposé (--ligament) : > 1 si la fissure est tortueuse
                               ou ramifiée
  énergie libérée              somme des W de breakHistory.txt (nJ)
  durée de propagation         de la première rupture à l'instant où la longueur rompue atteint 90 % de sa valeur
                               finale (ms)
  vitesse de propagation       0.9 x longueur finale / durée de propagation (m/s)

La longueur finale est celle de l'arrêt du calcul (5 ms après la chute de force) : si la fissure n'a pas fini de
traverser l'échantillon, elle est sous-estimée.

Produit : carte_ruptures.png, courbes_longueur_rompue.png, courbes_parametres_ruptures.png, ruptures.csv
"""

import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_maps import GRID, INK, INK2, SURFACE, heatmap, param_curves_figure, read_cases

SERIES = "#2a78d6"
LIGAMENT = 12.42e-3  # pointe effective de l'entaille en x = 1.181 mm, bord droit en x = 13.6 mm


def analyse_breaks(d):
    """Grandeurs de rupture d'un cas terminé (None sinon)."""
    status = os.path.join(d, "status.txt")
    if not os.path.exists(status) or not open(status).read().startswith("done"):
        return None
    ev = np.loadtxt(os.path.join(d, "breakEvol.txt"), comments="#", ndmin=2)
    bh = np.loadtxt(os.path.join(d, "breakHistory.txt"), comments="#", ndmin=2)
    if len(bh) == 0 or ev[-1, 2] <= 0.0:
        return None
    t, L = ev[:, 0], ev[:, 2]
    t1, Lf = bh[0, 0], L[-1]
    t90 = t[np.argmax(L >= 0.9 * Lf)]
    duree = t90 - t1
    return {
        "t": t, "L": L, "t_1re": t1, "t_fin": t[-1], "t_90": t90,
        "L_fin": Lf, "n_ruptures": len(bh), "E": float(bh[:, 9].sum()),
        "duree": duree, "vitesse": 0.9 * Lf / duree if duree > 0 else np.nan,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ligament", type=float, default=LIGAMENT,
                    help=f"longueur d'une fissure droite traversant le ligament, en m (défaut {LIGAMENT})")
    a = ap.parse_args()
    if not os.path.exists("cases.txt"):
        sys.exit("plot_breaks.py : cases.txt introuvable (lancer le script dans le répertoire des cas)")
    names, cases = read_cases()
    px, py = names[0], names[1]

    data = {}
    for d, p in cases:
        r = analyse_breaks(d)
        if r is not None:
            data[(p[px], p[py])] = (d, r)

    with open("ruptures.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cas", px, py, "liens_rompus", "longueur_rompue_mm", "longueur_sur_ligament", "energie_nJ",
                    "t_1re_rupture_ms", "duree_propagation_ms", "vitesse_m_par_s"])
        for d, p in cases:
            key = (p[px], p[py])
            if key not in data:
                w.writerow([d, p[px], p[py]] + [""] * 7)
                continue
            r = data[key][1]
            w.writerow([d, p[px], p[py], r["n_ruptures"], f"{r['L_fin'] * 1e3:.3f}", f"{r['L_fin'] / a.ligament:.3f}",
                        f"{r['E'] * 1e9:.2f}", f"{r['t_1re'] * 1e3:.3f}", f"{r['duree'] * 1e3:.3f}",
                        f"{r['vitesse']:.2f}"])

    xs = sorted({p[px] for _, p in cases})
    ys = sorted({p[py] for _, p in cases})
    plt.rcParams.update({"font.size": 10, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})

    # --- cartes
    grids = {k: np.full((len(ys), len(xs)), np.nan) for k in ("Lrel", "E", "duree", "vitesse")}
    for (x, y), (_, r) in data.items():
        j, i = ys.index(y), xs.index(x)
        grids["Lrel"][j, i] = r["L_fin"] / a.ligament
        grids["E"][j, i] = r["E"] * 1e9
        grids["duree"][j, i] = r["duree"] * 1e3
        grids["vitesse"][j, i] = r["vitesse"]
    cell = 0.55 * max(len(xs), 5) / 5
    fig, axs = plt.subplots(2, 2, figsize=(2 * 6.2 * cell, 2 * 5.2 * cell))
    heatmap(axs[0][0], xs, ys, grids["Lrel"], "Longueur rompue / ligament", "longueur relative", "{:.2f}")
    heatmap(axs[0][1], xs, ys, grids["E"], "Énergie élastique libérée par les ruptures", "E (nJ)", "{:.0f}")
    heatmap(axs[1][0], xs, ys, grids["duree"], "Durée de propagation (1re rupture → 90 %)", "durée (ms)", "{:.1f}")
    heatmap(axs[1][1], xs, ys, grids["vitesse"], "Vitesse moyenne de propagation", "v (m/s)", "{:.1f}")
    fs = 7 if len(xs) > 6 else 9
    for ax in axs.ravel():
        for txt in ax.texts:
            txt.set_fontsize(fs)
    fig.tight_layout()
    fig.savefig("carte_ruptures.png", dpi=150)

    param_curves_figure([(grids["Lrel"], "longueur rompue / ligament"), (grids["E"], "énergie libérée (nJ)"),
                         (grids["duree"], "durée de propagation (ms)"), (grids["vitesse"], "vitesse (m/s)")],
                        xs, ys, "courbes_parametres_ruptures.png", "Grandeurs de rupture en fonction de chaque paramètre")

    # --- courbes : longueur rompue en fonction du temps depuis la première rupture
    t_max = 1.1 * max((r["t_fin"] - r["t_1re"]) for _, r in data.values()) * 1e3
    l_max = 1.1 * max(max(r["L_fin"] for _, r in data.values()), a.ligament) * 1e3
    ys_desc = sorted(ys, reverse=True)
    fig, axs = plt.subplots(len(ys_desc), len(xs), figsize=(2.0 * len(xs) + 1.2, 1.6 * len(ys_desc) + 1.0),
                            sharex=True, sharey=True, squeeze=False)
    for j, y in enumerate(ys_desc):
        for i, x in enumerate(xs):
            ax = axs[j][i]
            ax.set_xlim(0, t_max)
            ax.set_ylim(0, l_max)
            ax.grid(True, color=GRID, linewidth=0.6)
            ax.tick_params(colors=INK2, labelsize=7, length=2)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color("#c9c8c0")
            ax.axhline(a.ligament * 1e3, color=INK2, linewidth=0.8, linestyle=":")
            if (x, y) not in data:
                ax.text(0.5, 0.5, "—", transform=ax.transAxes, ha="center", va="center", color=INK2)
                continue
            r = data[(x, y)][1]
            m = r["t"] >= r["t_1re"]
            ax.plot((r["t"][m] - r["t_1re"]) * 1e3, r["L"][m] * 1e3, color=SERIES, linewidth=1.5)
    for i, x in enumerate(xs):
        axs[-1][i].set_xlabel("t − t$_{1re}$ (ms)", color=INK2, fontsize=8)
        axs[0][i].set_title(f"$\\xi$ = {x:g}", color=INK, fontsize=10)
    for j, y in enumerate(ys_desc):
        axs[j][0].set_ylabel(f"$\\zeta$ = {y:g}\n\nL (mm)", color=INK, fontsize=9)
    handles = [plt.Line2D([], [], color=SERIES, linewidth=1.5, label="longueur rompue cumulée"),
               plt.Line2D([], [], color=INK2, linewidth=0.8, linestyle=":",
                          label=f"fissure droite à travers le ligament ({a.ligament * 1e3:.2f} mm)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Longueur rompue cumulée en fonction du temps depuis la première rupture", color=INK, fontsize=12)
    fig.tight_layout(rect=(0, 0.03, 1, 0.98))
    fig.savefig("courbes_longueur_rompue.png", dpi=150)

    print(f"{len(data)} cas exploités sur {len(cases)} -> carte_ruptures.png, courbes_longueur_rompue.png, "
          "courbes_parametres_ruptures.png, ruptures.csv")


if __name__ == "__main__":
    main()
