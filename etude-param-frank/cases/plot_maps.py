#!/usr/bin/env python3
"""Cartes de l'étude paramétrique dans l'espace (xi, zeta) : raideur et forces à la rupture.

À lancer dans le répertoire des cas (celui qui contient cases.txt) :

    python3 plot_maps.py

Pour chaque cas terminé (status.txt = done), on lit top.txt et bottom.txt (écrits par captureNodes toutes les
TCAPT secondes ; la période est lue dans l'input du cas : define TCAPT, ou à défaut define TOUT) :
  colonnes : x y (position moyenne des noeuds du mors)  Fx Fy (force totale)  fx fy (force moyenne)

  ouverture  D = (y_haut - y_bas) - (y_haut - y_bas)(t = 0)
  force      F = (Fy_bas - Fy_haut) / 2       (traction moyenne sur les deux mors)

  raideur K           pente de F(D) ajustée sur les points avant le pic où F < 50 % du pic (régime élastique)
  force au pic        maximum de F (échantillonné toutes les TCAPT secondes)
  force à la 1re rupture  force maximale atteinte jusqu'au premier lien rompu (instant lu dans breakHistory.txt)
                          -- égale au pic si la première rupture provoque la ruine (comportement fragile),
                          inférieure si des liens rompent bien avant le pic (endommagement progressif)

Produit : resultats.csv, carte_raideur.png, carte_forces_rupture.png.
"""

import csv
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize

# palette séquentielle à une teinte (bleu), clair -> foncé
CMAP = LinearSegmentedColormap.from_list("bleu", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])
INK, INK2, GRID, SURFACE = "#1f1e1b", "#5f5e57", "#e6e5df", "#fcfcfb"


def read_cases(path="cases.txt"):
    """Liste (répertoire, {paramètre: valeur}) lue dans cases.txt."""
    names, cases = None, []
    for line in open(path):
        if line.startswith("#"):
            tok = line[1:].split()
            if tok[:2] == ["id", "répertoire"]:
                names = tok[2:]
            continue
        tok = line.split()
        if len(tok) >= 2:
            cases.append((tok[1], dict(zip(names, map(float, tok[2:])))))
    return names, cases


def read_define(input_file, name):
    """Valeur numérique d'une ligne « define NAME valeur » (None si absente ou non numérique)."""
    for line in open(input_file):
        tok = line.split("#")[0].split()
        if len(tok) >= 3 and tok[0] == "define" and tok[1] == name:
            try:
                return float(tok[2])
            except ValueError:
                return None
    return None


def load_curve(d):
    """Courbe force-ouverture (D, F) d'un cas terminé, None sinon."""
    status = os.path.join(d, "status.txt")
    if not os.path.exists(status) or not open(status).read().startswith("done"):
        return None
    top, bot = np.loadtxt(os.path.join(d, "top.txt"), ndmin=2), np.loadtxt(os.path.join(d, "bottom.txt"), ndmin=2)
    n = min(len(top), len(bot))
    top, bot = top[:n], bot[:n]
    D = (top[:, 1] - bot[:, 1]) - (top[0, 1] - bot[0, 1])
    F = (bot[:, 3] - top[:, 3]) / 2.0
    return D, F


def capture_period(d):
    """Période d'écriture de top.txt / bottom.txt (define TCAPT, ou TOUT pour les études sans nstepPeriodCapture)."""
    inp = os.path.join(d, "input.txt")
    return read_define(inp, "TCAPT") or read_define(inp, "TOUT")


def analyse(d):
    """Grandeurs mesurées pour un cas (None si le cas n'est pas terminé ou pas exploitable)."""
    curve = load_curve(d)
    if curve is None:
        return None
    D, F = curve
    n = len(F)
    k = int(np.argmax(F))
    sel = (np.arange(n) <= k) & (F < 0.5 * F[k]) & (D > 0)
    if sel.sum() < 2:
        return None
    K = np.polyfit(D[sel], F[sel], 1)[0]

    r = {"K": K, "F_pic": F[k], "D_pic": D[k], "n_fit": int(sel.sum())}

    # force à la première rupture : maximum de F sur les sorties antérieures au premier lien rompu (une
    # interpolation à l'instant de rupture tomberait, dans les cas fragiles, au milieu de la chute de force)
    tout = capture_period(d)
    r["t_1re"] = r["F_1re"] = r["D_1re"] = math.nan
    hist = os.path.join(d, "breakHistory.txt")
    if tout and os.path.exists(hist):
        bh = np.loadtxt(hist, comments="#", ndmin=2)
        if len(bh):
            t1 = bh[0, 0]
            r["t_1re"] = t1
            before = np.arange(n) * tout <= t1
            r["F_1re"] = float(F[before].max()) if before.any() else math.nan
            r["D_1re"] = float(np.interp(t1, np.arange(n) * tout, D))

    # pic lissé relevé par l-hyphen (événement stopAfterStressDrop), pour comparaison
    r["F_pic_log"] = math.nan
    log = os.path.join(d, "log.txt")
    if os.path.exists(log):
        m = re.search(r"peak = ([0-9.eE+-]+)", open(log, errors="replace").read())
        if m:
            r["F_pic_log"] = float(m.group(1))
    return r


def heatmap(ax, xs, ys, Z, title, unit, fmt, norm=None):
    """Carte xi (abscisse) x zeta (ordonnée), cases régulières, valeur écrite dans chaque case."""
    if norm is None:
        finite = Z[np.isfinite(Z)]
        norm = Normalize(finite.min(), finite.max())
    masked = np.ma.masked_invalid(Z)
    cmap = CMAP.copy()
    cmap.set_bad("#eeede8")
    im = ax.imshow(masked, origin="lower", cmap=cmap, norm=norm, aspect="equal")
    for j in range(len(ys)):
        for i in range(len(xs)):
            v = Z[j, i]
            if np.isfinite(v):
                dark = norm(v) > 0.55
                ax.text(i, j, fmt.format(v), ha="center", va="center", fontsize=9,
                        color="#ffffff" if dark else INK)
            else:
                ax.text(i, j, "—", ha="center", va="center", fontsize=9, color=INK2)
    ax.set_xticks(range(len(xs)), [f"{x:g}" for x in xs])
    ax.set_yticks(range(len(ys)), [f"{y:g}" for y in ys])
    ax.set_xlabel(r"$\xi = k_b / (k_s\,\ell^2)$", color=INK)
    ax.set_ylabel(r"$\zeta = k_c / K_{cell}$", color=INK)
    ax.set_title(title, color=INK, fontsize=11)
    ax.tick_params(colors=INK2, length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(unit, color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=INK2)


def main():
    if not os.path.exists("cases.txt"):
        sys.exit("plot_maps.py : cases.txt introuvable (lancer le script dans le répertoire des cas)")
    names, cases = read_cases()
    if names is None or len(names) < 2:
        sys.exit("plot_maps.py : cases.txt doit avoir deux paramètres (xi, zeta)")
    px, py = names[0], names[1]  # XI en abscisse, ZETA en ordonnée

    results = []
    for d, p in cases:
        r = analyse(d)
        results.append((d, p, r))
        if r is None:
            print(f"  {d} : non terminé ou inexploitable, ignoré")

    # tableau
    with open("resultats.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cas", px, py, "K_N_par_m", "F_pic_mN", "D_pic_um", "F_1re_rupture_mN", "t_1re_rupture_ms",
                    "F_pic_log_mN", "points_ajustement"])
        for d, p, r in results:
            if r is None:
                w.writerow([d, p[px], p[py]] + [""] * 7)
            else:
                w.writerow([d, p[px], p[py], f"{r['K']:.1f}", f"{r['F_pic']*1e3:.3f}", f"{r['D_pic']*1e6:.2f}",
                            f"{r['F_1re']*1e3:.3f}", f"{r['t_1re']*1e3:.2f}", f"{r['F_pic_log']*1e3:.3f}",
                            r["n_fit"]])

    xs = sorted({p[px] for _, p, _ in results})
    ys = sorted({p[py] for _, p, _ in results})
    grids = {key: np.full((len(ys), len(xs)), np.nan) for key in ("K", "F_pic", "F_1re")}
    for d, p, r in results:
        if r is not None:
            j, i = ys.index(p[py]), xs.index(p[px])
            grids["K"][j, i] = r["K"]
            grids["F_pic"][j, i] = r["F_pic"] * 1e3
            grids["F_1re"][j, i] = r["F_1re"] * 1e3

    plt.rcParams.update({"font.size": 10, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    heatmap(ax, xs, ys, grids["K"], "Raideur de l'échantillon (pente élastique F – ouverture)", "K (N/m)",
            "{:.0f}")
    fig.tight_layout()
    fig.savefig("carte_raideur.png", dpi=150)

    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    both = np.concatenate([grids["F_pic"].ravel(), grids["F_1re"].ravel()])
    both = both[np.isfinite(both)]
    norm = Normalize(both.min(), both.max())  # même échelle de couleur pour les deux cartes de force
    heatmap(axs[0], xs, ys, grids["F_pic"], "Force au pic", "F (mN)", "{:.1f}", norm)
    heatmap(axs[1], xs, ys, grids["F_1re"], "Force maximale avant la première rupture de lien", "F (mN)", "{:.1f}",
            norm)
    fig.tight_layout()
    fig.savefig("carte_forces_rupture.png", dpi=150)

    ok = [r for _, _, r in results if r is not None]
    print(f"{len(ok)} cas exploités sur {len(results)} -> resultats.csv, carte_raideur.png, carte_forces_rupture.png")
    ecart = [abs(r["F_pic"] - r["F_pic_log"]) / r["F_pic_log"] for r in ok if np.isfinite(r["F_pic_log"])]
    if ecart:
        print(f"écart entre le pic échantillonné et le pic lissé de log.txt : {100*np.median(ecart):.1f} % "
              f"(médiane), {100*max(ecart):.1f} % (max) -- diminuer TCAPT dans input_model.txt pour affiner")


if __name__ == "__main__":
    main()
