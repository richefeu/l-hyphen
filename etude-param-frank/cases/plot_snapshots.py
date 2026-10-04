#!/usr/bin/env python3
"""Planches d'images see2 de tous les cas, disposées en carte dans l'espace (xi, zeta).

À lancer dans le répertoire des cas (celui qui contient cases.txt) :

    python3 plot_snapshots.py                     # toutes les vues
    python3 plot_snapshots.py --vue rupture       # une seule vue
    python3 plot_snapshots.py --see2 /chemin/see2 --taille 500 --jobs 4

Vues (une planche par vue, même disposition que plot_curves.py : xi en colonne, zeta croissant vers le haut) :
  rupture      faciès de rupture : dernier conf, liens rompus en rouge (see2 : show_crack_path)
  deformation  déformation eps_yy au pic de force, par rapport à conf0 (see2 : show_strain=eps_yy)
  contrainte   contrainte sig_yy au pic de force (see2 : show_stress=sig_yy, moyenne de Love-Weber, N/m)
  dirs_deformation  directions principales de déformation à l'initiation de la fissure (see2 : show_strain_dirs)
  dirs_contrainte   directions principales de contrainte à l'initiation de la fissure (see2 : show_stress_dirs)

Les conf sont choisis d'après leur temps (lu dans chaque fichier conf) :
  « au pic de force »  dernier conf avant le pic de force (force lue dans top.txt / bottom.txt, toutes les TCAPT
                       secondes) ;
  « à l'initiation »   conf sauvegardé à la première rupture de lien (événement saveConfAtBrokenLength 0 de
                       l'input) ; à défaut, dernier conf avant le premier lien rompu (breakHistory.txt).

Directions principales : un trait par direction et par cellule, épais pour la direction majeure, fin pour la
mineure, rouge en traction, bleu en compression, de longueur proportionnelle à la valeur principale. La même
valeur de référence (la plus grande valeur principale de tous les cas) est utilisée pour tous les cas : les
longueurs sont comparables d'un cas à l'autre. Avec --echelle par_cas, chaque cas est normalisé par sa propre
plus grande valeur principale (et chaque champ coloré par sa propre borne) : les orientations sont lisibles
même dans les cas peu chargés, mais les longueurs et les couleurs ne se comparent plus d'un cas à l'autre.

Les champs sont tracés avec la même échelle de couleur pour tous les cas : un premier passage de see2 relève la
borne automatique de chaque cas, puis toutes les images sont refaites avec la plus grande, et une seule barre de
couleur est ajoutée à la planche. Les images individuelles sont gardées dans chaque cas (see2_<vue>.png).

Produit : images_rupture.png, images_deformation.png, images_contrainte.png, images_dirs_deformation.png,
          images_dirs_contrainte.png (suffixe _par_cas avec --echelle par_cas)
"""

import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize

from plot_maps import INK, INK2, SURFACE, analyse, capture_period, load_curve, read_cases

# vues : réglages see2 et conf à rendre ("last" ou "pic")
VIEWS = {
    "rupture": {"conf": "last", "set": ["show_crack_path=1"], "titre": "Faciès de rupture (dernier conf, liens rompus en rouge)"},
    "deformation": {"conf": "pic", "set": ["show_strain=eps_yy", "refConf=0"], "champ": "strain",
                    "titre": r"Déformation $\varepsilon_{yy}$ au pic de force", "unite": r"$\varepsilon_{yy}$"},
    "contrainte": {"conf": "pic", "set": ["show_stress=sig_yy"], "champ": "stress",
                   "titre": r"Contrainte $\sigma_{yy}$ au pic de force", "unite": r"$\sigma_{yy}$ (N/m)"},
    "dirs_deformation": {"conf": "initiation", "set": ["show_strain_dirs=1", "refConf=0", "eScale=1"],
                         "dirs": "strain", "titre": "Directions principales de déformation à l'initiation de la fissure",
                         "unite": "déformation principale"},
    "dirs_contrainte": {"conf": "initiation", "set": ["show_stress_dirs=1", "eScale=1"], "dirs": "stress",
                        "titre": "Directions principales de contrainte à l'initiation de la fissure",
                        "unite": "contrainte principale (N/m)"},
}
COMMON = ["show_hud=0", "show_colorbar=0", "show_background=0", "show_cells=1", "show_contours=1"]

# palettes de see2 (TensorSphTable : bleu-blanc-rouge pour un champ signé ; TensorDevTable : blanc-jaune-rouge)
CMAP_SIGNED = LinearSegmentedColormap.from_list("see2_signe", [(40 / 255, 60 / 255, 200 / 255), (1, 1, 1),
                                                                (200 / 255, 30 / 255, 30 / 255)])
CMAP_POSITIVE = LinearSegmentedColormap.from_list("see2_positif", [(1, 1, 1), (1, 210 / 255, 60 / 255),
                                                                   (200 / 255, 30 / 255, 30 / 255)])


def conf_times(d):
    """{numéro: temps} des fichiers confN d'un cas (ligne « t ... » en tête de chaque conf)."""
    times = {}
    for f in os.listdir(d):
        m = re.fullmatch(r"conf(\d+)", f)
        if not m:
            continue
        with open(os.path.join(d, f)) as fh:
            for _ in range(40):
                tok = fh.readline().split()
                if len(tok) == 2 and tok[0] == "t":
                    times[int(m.group(1))] = float(tok[1])
                    break
    return times


def conf_before(times, t, tol=1e-12):
    """Numéro du dernier conf de temps <= t (None s'il n'y en a pas)."""
    ok = [(tt, n) for n, tt in times.items() if tt <= t + tol]
    return max(ok)[1] if ok else None


def run_see2(see2, case_dir, conf, sets, size, out):
    """Rend un conf avec see2 ; retourne (borne du champ, champ signé ?, temps) lus dans sa sortie."""
    cmd = [see2, str(conf), "--snapshot", out, "--size", f"{size}x{size}"]
    for s in COMMON + sets:
        cmd += ["--set", s]
    res = subprocess.run(cmd, cwd=case_dir, capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(os.path.join(case_dir, out)):
        raise RuntimeError(f"see2 a échoué dans {case_dir} :\n{res.stderr or res.stdout[-500:]}")
    bound, divergent, t = None, True, None
    m = re.search(r"field \S+ bound (\S+) divergent (\d)", res.stdout)
    if m:
        bound, divergent = float(m.group(1)), m.group(2) == "1"
    m = re.search(r"dirs \S+ max (\S+)", res.stdout)
    if m:
        bound = float(m.group(1))
    m = re.search(r"t = (\S+),", res.stdout)
    if m:
        t = float(m.group(1))
    return bound, divergent, t


def crop_box(images):
    """Rectangle commun contenant tous les pixels non blancs (les cas ont la même géométrie)."""
    r0, r1, c0, c1 = None, None, None, None
    for img in images:
        mask = np.any(img[:, :, :3] < 0.98, axis=2)
        rows, cols = np.where(mask.any(axis=1))[0], np.where(mask.any(axis=0))[0]
        if len(rows) == 0:
            continue
        r0 = rows[0] if r0 is None else min(r0, rows[0])
        r1 = rows[-1] if r1 is None else max(r1, rows[-1])
        c0 = cols[0] if c0 is None else min(c0, cols[0])
        c1 = cols[-1] if c1 is None else max(c1, cols[-1])
    pad = 4
    h, w = images[0].shape[:2]
    return max(r0 - pad, 0), min(r1 + pad + 1, h), max(c0 - pad, 0), min(c1 + pad + 1, w)


def make_view(name, cases, px, py, see2, size, jobs, common=True):
    view = VIEWS[name]
    out = f"see2_{name}.png"

    # conf à rendre pour chaque cas
    todo = []
    for d, p in cases:
        curve = load_curve(d)
        if curve is None:
            print(f"  {d} : non terminé, case vide")
            continue
        times = conf_times(d)
        if not times:
            print(f"  {d} : aucun conf, case vide")
            continue
        if view["conf"] == "last":
            conf = max(times)
        elif view["conf"] == "pic":
            tcapt = capture_period(d)
            conf = conf_before(times, int(np.argmax(curve[1])) * tcapt) if tcapt else None
        else:  # initiation : conf sauvegardé à la première rupture, sinon dernier conf avant celle-ci
            r = analyse(d)
            if r is None or not np.isfinite(r["t_1re"]):
                print(f"  {d} : pas de lien rompu, case vide")
                continue
            t1 = r["t_1re"]
            after = [(tt, n) for n, tt in times.items() if t1 - 1e-12 <= tt <= t1 + 2e-5]
            conf = min(after)[1] if after else conf_before(times, t1)
        if conf is None:
            print(f"  {d} : pas de conf à l'instant voulu, case vide")
            continue
        todo.append((d, p, conf))

    def render(item, extra):
        d, p, conf = item
        return run_see2(see2, d, conf, view["set"] + extra, size, out)

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        first = list(pool.map(lambda it: render(it, []), todo))
        cmax, divergent, dmax = None, True, None
        if "dirs" in view and common:
            # même longueur de référence pour tous les cas : la plus grande valeur principale
            dmax = max(b for b, _, _ in first if b)
            key = "strain_dirsMax" if view["dirs"] == "strain" else "stress_dirsMax"
            first = list(pool.map(lambda it: render(it, [f"{key}={dmax}"]), todo))
        if "champ" in view and common:
            # même échelle pour tous les cas : on refait les images avec la plus grande borne automatique
            bounds = [b for b, _, _ in first if b]
            cmax, divergent = max(bounds), first[0][1]
            key = "strain_colorMax" if view["champ"] == "strain" else "stress_colorMax"
            first = list(pool.map(lambda it: render(it, [f"{key}={cmax}"]), todo))

    images = {(p[px], p[py]): (plt.imread(os.path.join(d, out)), t) for (d, p, _), (_, _, t) in zip(todo, first)}
    r0, r1, c0, c1 = crop_box([img for img, _ in images.values()])

    xs = sorted({p[px] for _, p in cases})
    ys = sorted({p[py] for _, p in cases}, reverse=True)  # zeta croissant vers le haut
    plt.rcParams.update({"font.size": 9, "figure.facecolor": SURFACE})
    fig, axs = plt.subplots(len(ys), len(xs), figsize=(2.3 * len(xs) + 1.6, 2.3 * len(ys) + 1.2), squeeze=False)
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            ax = axs[j][i]
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            if (x, y) not in images:
                ax.text(0.5, 0.5, "—", transform=ax.transAxes, ha="center", va="center", color=INK2)
                continue
            img, t = images[(x, y)]
            ax.imshow(img[r0:r1, c0:c1])
            if t is not None:
                ax.set_xlabel(f"t = {t * 1e3:.1f} ms", fontsize=7, color=INK2, labelpad=2)
    for i, x in enumerate(xs):
        axs[0][i].set_title(f"$\\xi$ = {x:g}", color=INK, fontsize=10)
    for j, y in enumerate(ys):
        axs[j][0].set_ylabel(f"$\\zeta$ = {y:g}", color=INK, fontsize=10)

    fig.suptitle(view["titre"] + ("" if common or not ("champ" in view or "dirs" in view)
                                  else " — échelle propre à chaque cas"), color=INK, fontsize=12)
    # mise en page des cases d'abord, barre de couleur ensuite (dans la marge réservée à droite)
    fig.tight_layout(rect=(0, 0, 0.9 if cmax is not None else 1.0, 0.97))
    if cmax is not None:
        cax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
        norm = Normalize(-cmax if divergent else 0.0, cmax)
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP_SIGNED if divergent else CMAP_POSITIVE), cax=cax)
        cb.set_label(view["unite"] + "  (échelle commune)", color=INK2)
        cb.outline.set_visible(False)
        cb.ax.tick_params(colors=INK2, labelsize=8)
    if "dirs" in view:
        handles = [
            plt.Line2D([], [], color="#d90d0d", linewidth=3.0, label="majeure, traction"),
            plt.Line2D([], [], color="#d90d0d", linewidth=1.0, label="mineure, traction"),
            plt.Line2D([], [], color="#0d26e6", linewidth=3.0, label="majeure, compression"),
            plt.Line2D([], [], color="#0d26e6", linewidth=1.0, label="mineure, compression"),
        ]
        title = (f"longueur ∝ |{view['unite']}| ; demi-trait d'un rayon de cellule pour {dmax:.3g} "
                 f"(plus grande valeur de tous les cas)" if dmax is not None else
                 "longueur ∝ valeur principale, normalisée dans chaque cas par sa plus grande valeur")
        fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9, title=title,
                   title_fontsize=8)
        fig.tight_layout(rect=(0, 0.07, 1.0, 0.97))
    suffix = "" if common or not ("champ" in view or "dirs" in view) else "_par_cas"
    fig.savefig(f"images_{name}{suffix}.png", dpi=150)
    scale = f", échelle commune ±{cmax:.4g}" if cmax else (f", longueur de référence {dmax:.4g}" if dmax else "")
    print(f"images_{name}{suffix}.png : {len(images)} cas{scale}")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vue", choices=list(VIEWS) + ["toutes"], default="toutes")
    ap.add_argument("--see2", default=os.path.normpath(os.path.join(here, "..", "..", "see2")),
                    help="exécutable see2 (défaut : see2 à la racine du dépôt)")
    ap.add_argument("--taille", type=int, default=500, help="taille des images see2 en pixels (défaut 500)")
    ap.add_argument("--jobs", type=int, default=4, help="nombre de see2 lancés en parallèle (défaut 4)")
    ap.add_argument("--echelle", choices=["commune", "par_cas"], default="commune",
                    help="échelle des champs et des directions : commune à tous les cas (défaut) ou propre à chacun")
    a = ap.parse_args()

    if not os.path.exists("cases.txt"):
        sys.exit("plot_snapshots.py : cases.txt introuvable (lancer le script dans le répertoire des cas)")
    if not os.access(a.see2, os.X_OK):
        sys.exit(f"plot_snapshots.py : see2 introuvable ({a.see2}) ; le compiler à la racine (make see2) ou --see2")
    see2 = os.path.abspath(a.see2)
    names, cases = read_cases()
    px, py = names[0], names[1]
    for name in (VIEWS if a.vue == "toutes" else [a.vue]):
        make_view(name, cases, px, py, see2, a.taille, a.jobs, a.echelle == "commune")


if __name__ == "__main__":
    main()
