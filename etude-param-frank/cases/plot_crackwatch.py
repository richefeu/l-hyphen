#!/usr/bin/env python3
"""Fissure principale et ramifications (crackwatch) dans l'espace (xi, zeta).

À lancer dans le répertoire des cas (celui qui contient cases.txt) :

    python3 plot_crackwatch.py [--crackwatch CHEMIN]

Pour chaque cas terminé, crackwatch (racine du dépôt) est lancé dans le répertoire du cas si ses sorties manquent ou
sont plus anciennes que breakHistory.txt. Il construit le graphe des liens rompus (positions lues dans conf0), le
découpe en fissures (composantes connexes) et, pour chacune, détermine son chemin principal : le plus long des plus
courts chemins entre deux extrémités. Ce chemin peut relier deux pointes de ramifications sans passer par l'entaille ;
le script lit donc le graphe complet (crackwatch_graph.txt : une ligne par segment rompu, avec l'instant de sa rupture)
et définit la fissure principale comme le plus court chemin de la pointe de l'entaille (--pointe) à l'extrémité de la
même fissure la plus éloignée de cette pointe. Le reste des liens rompus forme les ramifications (et les fissures
secondaires).

Grandeurs :
  tortuosité            longueur de la fissure principale / distance entre ses extrémités (1 = fissure droite)
  part des ramifications  1 - longueur de la fissure principale / longueur rompue totale (toutes fissures)
  vitesse               pente de la longueur rompue de la fissure principale en fonction du temps, ajustée entre 20 % et
                        60 % de sa longueur finale (m/s) : avant les arrêts observés vers 80 %
  délai d'amorçage      de la première rupture à l'instant où la fissure principale atteint 20 % de sa longueur (ms)

Produit : carte_crackwatch.png, courbes_parametres_crackwatch.png, courbes_fissure_principale.png, facies_crackwatch.png,
          crackwatch.csv
"""

import argparse
import csv
import heapq
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_maps import GRID, INK, INK2, SURFACE, heatmap, param_curves_figure, read_cases

SERIES = "#2a78d6"
MAIN = "#0d366b"    # fissure principale
BRANCH = "#e0812f"  # ramifications et fissures secondaires
BOX = (13.6e-3, 16e-3)  # dimensions de l'échantillon (m)
TIP = (1.181e-3, 8e-3)  # pointe effective de l'entaille de echantillon/nodefile.txt


def run_crackwatch(crackwatch, d):
    """Lance crackwatch dans le cas d si ses sorties manquent ou sont plus anciennes que breakHistory.txt."""
    hist = os.path.join(d, "breakHistory.txt")
    out = os.path.join(d, "crackwatch_graph.txt")
    if os.path.exists(out) and os.path.getmtime(out) >= os.path.getmtime(hist):
        return True
    res = subprocess.run([crackwatch], cwd=d, capture_output=True, text=True)
    return res.returncode == 0 and os.path.exists(out)


def shortest_paths(adj, src):
    """Dijkstra : distances et prédécesseurs depuis src (adj[n] = [(voisin, longueur, arête)])."""
    dist, prev, heap = {src: 0.0}, {src: None}, [(0.0, src)]
    while heap:
        dn, n = heapq.heappop(heap)
        if dn > dist[n]:
            continue
        for m, w, e in adj[n]:
            if dn + w < dist.get(m, np.inf):
                dist[m], prev[m] = dn + w, (n, e)
                heapq.heappush(heap, (dn + w, m))
    return dist, prev


def analyse_crack(d, crackwatch, tip):
    """Grandeurs de la fissure principale d'un cas terminé (None sinon).

    Fissure principale : chemin le plus court, dans le graphe des liens rompus, de la pointe de l'entaille (nœud le
    plus proche de tip) à l'extrémité de la même fissure la plus éloignée (à vol d'oiseau) de cette pointe.
    """
    status = os.path.join(d, "status.txt")
    if not os.path.exists(status) or not open(status).read().startswith("done"):
        return None
    if not os.path.exists(os.path.join(d, "breakHistory.txt")) or not run_crackwatch(crackwatch, d):
        return None
    g = np.loadtxt(os.path.join(d, "crackwatch_graph.txt"), comments="#", ndmin=2)
    if len(g) == 0:
        return None
    # nœuds globaux (fissure, nœud) -> position ; arêtes -> (instant, longueur)
    pos, adj = {}, {}
    for e, (c, a, b, t, w, xa, ya, xb, yb) in enumerate(g):
        na, nb = (int(c), int(a)), (int(c), int(b))
        pos[na], pos[nb] = (xa, ya), (xb, yb)
        adj.setdefault(na, []).append((nb, w, e))
        adj.setdefault(nb, []).append((na, w, e))
    keys = list(pos)
    P = np.array([pos[k] for k in keys])
    src = keys[int(np.argmin(np.hypot(P[:, 0] - tip[0], P[:, 1] - tip[1])))]
    dist, prev = shortest_paths(adj, src)
    ends = [n for n in dist if len(adj[n]) == 1 and n != src] or [n for n in dist if n != src]
    if not ends:
        return None
    far = max(ends, key=lambda n: np.hypot(pos[n][0] - pos[src][0], pos[n][1] - pos[src][1]))
    edges, n = [], far
    while prev[n] is not None:
        n, e = prev[n]
        edges.append(e)
    main = g[edges]
    L_main, L_tot = main[:, 4].sum(), g[:, 4].sum()
    chord = np.hypot(pos[far][0] - pos[src][0], pos[far][1] - pos[src][1])
    t1 = g[:, 3].min()  # première rupture (toutes fissures)
    comp = g[g[:, 0] == src[0]]
    nodes = {}
    for _, a, b, *_ in comp:
        nodes[a] = nodes.get(a, 0) + 1
        nodes[b] = nodes.get(b, 0) + 1

    # croissance de la fissure principale : longueur cumulée de ses segments, par instant de rupture
    order = np.argsort(main[:, 3])
    t = main[order, 3]
    Lc = np.cumsum(main[order, 4])
    sel = (Lc >= 0.2 * L_main) & (Lc <= 0.6 * L_main)
    v = np.polyfit(t[sel], Lc[sel], 1)[0] if sel.sum() >= 3 else np.nan
    t20 = t[np.argmax(Lc >= 0.2 * L_main)]
    return {
        "n_fissures": len(np.unique(g[:, 0])), "L_main": L_main, "L_tot": L_tot, "chord": chord,
        "tortuosite": L_main / chord if chord > 0 else np.nan, "ramif": 1.0 - L_main / L_tot,
        "bifurcations": sum(1 for k in nodes.values() if k >= 3), "extremites": sum(1 for k in nodes.values() if k == 1),
        "vitesse": v, "delai": t20 - t1, "t1": t1, "t": t, "Lc": Lc, "g": g, "main_edges": edges,
    }


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crackwatch", default=os.path.normpath(os.path.join(here, "..", "..", "crackwatch")),
                    help="exécutable crackwatch (défaut : crackwatch à la racine du dépôt)")
    ap.add_argument("--pointe", type=float, nargs=2, default=TIP, metavar=("X", "Y"),
                    help=f"position de la pointe de l'entaille, en m (défaut {TIP[0]} {TIP[1]})")
    a = ap.parse_args()
    if not os.path.exists("cases.txt"):
        sys.exit("plot_crackwatch.py : cases.txt introuvable (lancer le script dans le répertoire des cas)")
    if not os.access(a.crackwatch, os.X_OK):
        sys.exit(f"plot_crackwatch.py : crackwatch introuvable ({a.crackwatch}) ; le compiler à la racine (make)")
    crackwatch = os.path.abspath(a.crackwatch)
    names, cases = read_cases()
    px, py = names[0], names[1]

    data = {}
    for d, p in cases:
        r = analyse_crack(d, crackwatch, a.pointe)
        if r is not None:
            data[(p[px], p[py])] = (d, r)

    with open("crackwatch.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cas", px, py, "fissures", "longueur_principale_mm", "longueur_totale_mm", "corde_mm",
                    "tortuosite", "part_ramifications", "bifurcations", "extremites", "vitesse_m_par_s", "delai_ms"])
        for d, p in cases:
            key = (p[px], p[py])
            if key not in data:
                w.writerow([d, p[px], p[py]] + [""] * 10)
                continue
            r = data[key][1]
            w.writerow([d, p[px], p[py], r["n_fissures"], f"{r['L_main'] * 1e3:.3f}", f"{r['L_tot'] * 1e3:.3f}",
                        f"{r['chord'] * 1e3:.3f}", f"{r['tortuosite']:.3f}", f"{r['ramif']:.3f}", r["bifurcations"],
                        r["extremites"], f"{r['vitesse']:.2f}", f"{r['delai'] * 1e3:.3f}"])

    xs = sorted({p[px] for _, p in cases})
    ys = sorted({p[py] for _, p in cases})
    plt.rcParams.update({"font.size": 10, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})

    keys = ("tortuosite", "ramif", "vitesse", "delai")
    grids = {k: np.full((len(ys), len(xs)), np.nan) for k in keys}
    for (x, y), (_, r) in data.items():
        j, i = ys.index(y), xs.index(x)
        grids["tortuosite"][j, i] = r["tortuosite"]
        grids["ramif"][j, i] = r["ramif"]
        grids["vitesse"][j, i] = r["vitesse"]
        grids["delai"][j, i] = r["delai"] * 1e3

    # --- cartes
    cell = 0.55 * max(len(xs), 5) / 5
    fig, axs = plt.subplots(2, 2, figsize=(2 * 6.2 * cell, 2 * 5.2 * cell))
    heatmap(axs[0][0], xs, ys, grids["tortuosite"], "Tortuosité de la fissure principale", "chemin / corde", "{:.2f}")
    heatmap(axs[0][1], xs, ys, grids["ramif"], "Part des ramifications dans la longueur rompue", "fraction", "{:.2f}")
    heatmap(axs[1][0], xs, ys, grids["vitesse"], "Vitesse de la fissure principale (20 → 60 %)", "v (m/s)", "{:.1f}")
    heatmap(axs[1][1], xs, ys, grids["delai"], "Délai d'amorçage (1re rupture → 20 %)", "délai (ms)", "{:.2f}")
    fs = 7 if len(xs) > 6 else 9
    for ax in axs.ravel():
        for txt in ax.texts:
            txt.set_fontsize(fs)
    fig.tight_layout()
    fig.savefig("carte_crackwatch.png", dpi=150)
    plt.close(fig)

    # --- courbes en fonction de chaque paramètre
    param_curves_figure([(grids["tortuosite"], "tortuosité (chemin / corde)"),
                         (grids["ramif"], "part des ramifications"),
                         (grids["vitesse"], "vitesse fissure principale (m/s)"),
                         (grids["delai"], "délai d'amorçage (ms)")],
                        xs, ys, "courbes_parametres_crackwatch.png",
                        "Fissure principale et ramifications (crackwatch) en fonction de chaque paramètre")

    # --- longueur de la fissure principale en fonction du temps depuis la première rupture
    t_max = 1.1 * max((r["t"][-1] - r["t1"]) for _, r in data.values()) * 1e3
    l_max = 1.1 * max(r["L_main"] for _, r in data.values()) * 1e3
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
            if (x, y) not in data:
                ax.text(0.5, 0.5, "—", transform=ax.transAxes, ha="center", va="center", color=INK2)
                continue
            r = data[(x, y)][1]
            ax.step(np.r_[0.0, (r["t"] - r["t1"]) * 1e3], np.r_[0.0, r["Lc"] * 1e3], where="post", color=SERIES,
                    linewidth=1.5)
    for i, x in enumerate(xs):
        axs[-1][i].set_xlabel("t − t$_{1re}$ (ms)", color=INK2, fontsize=8)
        axs[0][i].set_title(f"$\\xi$ = {x:g}", color=INK, fontsize=10)
    for j, y in enumerate(ys_desc):
        axs[j][0].set_ylabel(f"$\\zeta$ = {y:g}\n\nL (mm)", color=INK, fontsize=9)
    fig.suptitle("Longueur rompue de la fissure principale (sans ramifications) en fonction du temps", color=INK,
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig("courbes_fissure_principale.png", dpi=150)
    plt.close(fig)

    # --- faciès : fissure principale et ramifications, dans la configuration initiale (conf0)
    # cadrage vertical commun : étendue des fissures issues de l'entaille (ramifications comprises), plus une marge ;
    # les petites fissures secondaires éloignées (près des mors) peuvent sortir du cadre
    ally = np.concatenate([r["g"][r["g"][:, 0] == r["g"][r["main_edges"][0], 0]][:, [6, 8]].ravel()
                           for _, r in data.values()]) * 1e3
    y0, y1 = max(ally.min() - 0.5, 0.0), min(ally.max() + 0.5, BOX[1] * 1e3)
    w = 1.5
    fig, axs = plt.subplots(len(ys_desc), len(xs), figsize=(w * len(xs) + 1.0,
                                                            w * (y1 - y0) / (BOX[0] * 1e3) * len(ys_desc) + 1.2),
                            squeeze=False)
    for j, y in enumerate(ys_desc):
        for i, x in enumerate(xs):
            ax = axs[j][i]
            ax.set_xlim(-0.2, BOX[0] * 1e3 + 0.2)
            ax.set_ylim(y0, y1)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for side in ax.spines.values():
                side.set_visible(False)
            ax.add_patch(plt.Rectangle((0, 0), BOX[0] * 1e3, BOX[1] * 1e3, fill=False, edgecolor="#c9c8c0",
                                       linewidth=0.8))
            ax.plot([0, a.pointe[0] * 1e3], [a.pointe[1] * 1e3] * 2, color=INK2, linewidth=1.2)
            if (x, y) not in data:
                ax.text(0.5, 0.5, "—", transform=ax.transAxes, ha="center", va="center", color=INK2)
                continue
            r = data[(x, y)][1]
            g = r["g"] * 1e3
            is_main = np.zeros(len(g), dtype=bool)
            is_main[r["main_edges"]] = True
            for m, color, lw in ((~is_main, BRANCH, 0.9), (is_main, MAIN, 1.3)):
                seg = g[m]
                X = np.column_stack([seg[:, 5], seg[:, 7], np.full(len(seg), np.nan)]).ravel()
                Y = np.column_stack([seg[:, 6], seg[:, 8], np.full(len(seg), np.nan)]).ravel()
                ax.plot(X, Y, color=color, linewidth=lw, solid_capstyle="round")
    for i, x in enumerate(xs):
        axs[0][i].set_title(f"$\\xi$ = {x:g}", color=INK, fontsize=10)
    for j, y in enumerate(ys_desc):
        axs[j][0].set_ylabel(f"$\\zeta$ = {y:g}", color=INK, fontsize=10)
    handles = [plt.Line2D([], [], color=MAIN, linewidth=1.5, label="fissure principale (depuis la pointe de l'entaille)"),
               plt.Line2D([], [], color=BRANCH, linewidth=1.5, label="ramifications et fissures secondaires"),
               plt.Line2D([], [], color=INK2, linewidth=1.2, label="entaille")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=10)
    fig.suptitle("Faciès de rupture (crackwatch) : fissure principale et ramifications, configuration initiale\n"
                 f"(largeur complète de l'échantillon, hauteur limitée à y = {y0:.1f} – {y1:.1f} mm)", color=INK, fontsize=12)
    fig.tight_layout(rect=(0, 0.025, 1, 0.98))
    fig.savefig("facies_crackwatch.png", dpi=150)
    plt.close(fig)

    print(f"{len(data)} cas exploités sur {len(cases)} -> carte_crackwatch.png, courbes_parametres_crackwatch.png, "
          f"courbes_fissure_principale.png, facies_crackwatch.png, crackwatch.csv")


if __name__ == "__main__":
    main()
