# Étude paramétrique (xi, zeta) sur des échantillons de type *Frank-ouverture-fissure*

Reprise complète de l'étude de `etude-param-model/` avec un échantillon de mêmes dimensions et de même nombre de
cellules que *Frank-ouverture-fissure* (13.6 mm × 16 mm, environ 2200 cellules), et une grille plus fine :

```
xi, zeta ∈ {0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 7.5, 10}      soit 11 × 11 = 121 cas
```

Les outils sont les mêmes que dans `etude-param-model/` (dont le `README.md` détaille chacun d'eux) ; ce
document donne la marche à suivre complète et ce qui change pour cette étude.

| Étape | Fichiers |
|---|---|
| 1. Échantillon | `echantillon/params.txt` → `nodegen` → `echantillon/nodefile.txt` |
| 2. Modèle d'input | `input_model.txt` (paramètres `XI`, `ZETA`) |
| 3. Génération des cas | `plan.txt` → `generate_cases` → `cases/` (121 répertoires) |
| 4. Lancement | `run_cases.sh` (parallèle, reprise) |
| Dépouillement | `cases/plot_maps.py`, `cases/plot_curves.py`, `cases/plot_snapshots.py`, `cases/plot_breaks.py`, `cases/plot_crackwatch.py` |

---

## 1. Échantillon

```bash
cd nodeGen && make && cd ../etude-param-frank/echantillon
../../nodeGen/nodegen params.txt
```

| | *Frank-ouverture-fissure* (`cellPrepro`) | cette étude (`nodegen`) |
|---|---|---|
| domaine | 13.6 mm × 16 mm | 13.6 mm × 16 mm |
| cellules / nœuds | 2155 / 12798 | 2222 / 13084 |
| diamètre équivalent des cellules | 0.353 mm | 0.353 mm (`cellSize`) |
| pré-entaille | 1.2 mm environ, cellules retirées | 1.2 mm, nette, depuis le milieu du bord gauche |
| `l_min / l_moy` | 0.012 | 0.40 (pas de barre courte) |
| cellules à 6 côtés (toutes les cellules) | 66 % | 88 % (94 % des cellules intérieures ; `lattice hex`, `disorder 0.3`) |

Le pavage est le même que celui du tutoriel de `nodegen` (`nodeGen/doc/tutorial_params.txt`), placé à
l'origine. La pointe effective de l'entaille est en (1.181 mm, 8 mm). Les mors attrapent la rangée de nœuds du
bord : 82 nœuds en bas, 79 en haut.

![Échantillon](echantillon/echantillon.png)

## 2. Modèle d'input `input_model.txt`

Même construction que dans `etude-param-model` : seules les lignes `define XI` et `define ZETA` changent d'un
cas à l'autre, et l-hyphen calcule

```
kb = xi ks l²        Kcell = 12 xi ks / (1 + 9 xi)        kc = zeta Kcell
```

avec `ks = 2000 N/m` et `l = 1.986276e-4 m` (longueur moyenne des barres de cet échantillon). `kb` sert à la
flexion des parois, `kc` aux liens cohésifs et au contact. Pour `xi = zeta = 10`, ce sont exactement les valeurs
de *Frank-ouverture-fissure* (à `l` près).

`dt = 5e-8 s` est stable sur toute la grille (`dt_crit/dt ≥ 16`, flexion pour `xi = 10`).

**Ce qui change par rapport à `etude-param-model` : les sorties.** Avec 2222 cellules, un conf pèse 3.3 Mo et un
`sample*.svg` environ 2 Mo : avec les réglages de la première étude, l'étude complète aurait occupé plusieurs
dizaines de Go. Le modèle utilise donc le mot-clé `nstepPeriodCapture` (ajouté à l-hyphen pour cette étude), qui
découple l'écriture de `top.txt` / `bottom.txt` de celle des svg :

| Sortie | Période | Rôle |
|---|---|---|
| `sample*.svg` | aucune (`nstepPeriodSVG 0`) | |
| `top.txt`, `bottom.txt` | `TCAPT = 0.1 ms` (`nstepPeriodCapture`) | courbes force – ouverture ; 5 fois plus fines que dans la première étude, d'où un pic mieux résolu |
| `conf*` | `TCONF = 1 ms`, plus un conf à la première rupture (`event saveConfAtBrokenLength 0`) | images see2 |

L'arrêt se fait 5 ms après une chute de 30 % de la force du mors haut (`event stopAfterStressDrop`), au lieu de
1 ms dans la première étude : la fissure a le temps de traverser l'échantillon (environ 3 ms dans *Frank*), ce
qui donne des faciès de rupture complets. `TFINAL = 0.15 s` n'est qu'une durée de sécurité.

Calibration sur les quatre coins et le centre de la grille (4 calculs simultanés sur cette machine) :

| xi | zeta | pic de force | chute détectée | arrêt | temps de calcul | disque |
|---|---|---|---|---|---|---|
| 0.1 | 0.1 | 8.3 mN | 50.7 ms | 55.7 ms | 19 min | 192 Mo |
| 0.1 | 10 | 44.7 mN | 59.7 ms | 64.7 ms | 22 min | 219 Mo |
| 1 | 1 | 24.0 mN | 23.5 ms | 28.5 ms | 11 min | 104 Mo |
| 10 | 0.1 | 12.1 mN | 35.6 ms | 40.6 ms | 14 min | 147 Mo |
| 10 | 10 | 42.5 mN | 28.4 ms | 33.4 ms | 11 min | 122 Mo |

Tous les cas s'arrêtent par l'événement, bien avant `TFINAL`. Pour l'étude complète, compter **10 à 25 min de
calcul par cas** et **environ 150 Mo par cas, soit 15 à 25 Go pour les 121 cas** (essentiellement les conf).
Sur une machine à 4 cœurs performance, l'étude prend de l'ordre de 8 h avec `-j 4`.

## 3. Générer les cas

```bash
cd etude-param-frank
make                          # compile generate_cases
./generate_cases              # crée cases/ : 121 répertoires XI<xi>_ZETA<zeta>/ (input.txt + nodefile.txt)
```

Un cas existant n'est pas modifié (ses résultats sont préservés) ; `--force` réécrit les `input.txt` sans toucher
aux résultats ; `--dry-run` affiche la liste des cas.

## 4. Lancer les calculs

Sur la machine de calcul (les cas sont ignorés par git, on les régénère sur place) :

```bash
git clone <dépôt> l-hyphen && cd l-hyphen
make                                     # compile run et see2
cd etude-param-frank
make && ./generate_cases
nohup ./run_cases.sh -j <nombre de cœurs> > run_cases.log 2>&1 &
./run_cases.sh -s                        # état de l'étude à tout moment
```

Sous macOS, empêcher la mise en veille pendant les calculs (sinon ils sont suspendus) :

```bash
nohup caffeinate -i ./run_cases.sh -j 4 > run_cases.log 2>&1 &
```

**Important : compiler l-hyphen à partir d'un répertoire propre** (`make clean && make` à la racine) si des
fichiers objets d'une version antérieure existent. Avant l'ajout des dépendances aux en-têtes dans le `Makefile`,
un objet non recompilé (`Event.o`) lisait les membres de `Lhyphen` à de mauvaises positions après l'ajout de
`nstepPeriodCapture` : les conf des événements écrasaient `conf0` et `conf1`.

Pour reprendre après un problème (interruption, redémarrage, échec), relancer la même commande : les cas
terminés sont sautés, les cas interrompus ou en échec sont relancés depuis le début (au plus 3 tentatives,
option `-m`). Voir `etude-param-model/README.md`, étape 4, pour le détail (`status.txt`, options, détection des
échecs).

## Dépouillement

```bash
cd etude-param-frank/cases
python3 plot_maps.py          # carte_raideur.png, carte_forces_rupture.png, courbes_parametres.png, resultats.csv
python3 plot_curves.py        # courbes_force_ouverture.png
python3 plot_snapshots.py     # images_rupture.png, images_deformation.png, images_contrainte.png,
                              # images_dirs_deformation.png, images_dirs_contrainte.png (see2 compilé requis)
python3 plot_breaks.py        # carte_ruptures.png, courbes_longueur_rompue.png, courbes_parametres_ruptures.png,
                              # ruptures.csv
python3 plot_crackwatch.py    # carte_crackwatch.png, courbes_parametres_crackwatch.png,
                              # courbes_fissure_principale.png, facies_crackwatch.png, crackwatch.csv
                              # (crackwatch compilé requis)
```

Chaque grandeur des cartes est aussi tracée en courbes (`courbes_parametres.png`,
`courbes_parametres_ruptures.png`) : en fonction de `xi`, une courbe par valeur de `zeta`, et en fonction de
`zeta`, une courbe par valeur de `xi`. La couleur des courbes suit la valeur du paramètre fixé (une teinte, du
clair au foncé, échelle logarithmique), lue sur la barre de couleur de chaque panneau.

`plot_breaks.py` exploite `breakHistory.txt` (une ligne par face collée rompue, avec l'énergie élastique `W` de
ses points de colle) et `breakEvol.txt` (longueur rompue cumulée toutes les 5 µs ; la 2e colonne de ce fichier est
la somme des taux de restitution `G` des faces rompues, en J/m, et non une énergie) :

| Carte | Définition |
|---|---|
| longueur rompue / ligament | longueur finale des faces rompues divisée par celle d'une fissure droite de la pointe de l'entaille au bord droit (12.42 mm, option `--ligament`) : > 1 pour une fissure tortueuse ou ramifiée |
| énergie libérée | somme des `W` de `breakHistory.txt` (nJ) |
| durée de propagation | de la première rupture à 90 % de la longueur rompue finale (ms) |
| vitesse de propagation | 0.9 × longueur finale / durée de propagation (m/s) |

La grille de courbes donne la longueur rompue cumulée en fonction du temps depuis la première rupture, avec les
mêmes axes pour tous les cas et, en pointillé, la longueur d'une fissure droite. Elle montre les cas où la fissure
s'arrête après une première rupture avant de repartir (par exemple `xi = 2`, `zeta ≥ 5` : la durée de propagation
inclut alors cet arrêt). La longueur finale est celle de l'arrêt du calcul (5 ms après la chute de force).

### Fissure principale et ramifications (`plot_crackwatch.py`)

Les grandeurs de `plot_breaks.py` mélangent la fissure principale et ses ramifications, dont l'apparition varie
beaucoup d'un cas à l'autre : elles sont bruitées. `plot_crackwatch.py` les sépare à l'aide de `crackwatch`
(racine du dépôt, `make crackwatch`), lancé dans chaque cas terminé si ses sorties manquent ou sont plus anciennes
que `breakHistory.txt`. `crackwatch` construit le graphe des faces rompues (positions lues dans `conf0`) et écrit
`crackwatch_summary.txt` (une ligne par fissure), `crackwatch_main.txt` (son propre chemin principal) et
`crackwatch_graph.txt` (tous les segments rompus, avec l'instant de leur rupture).

Le chemin principal de `crackwatch` (le plus long des plus courts chemins entre deux extrémités) peut relier deux
pointes de ramifications sans passer par l'entaille (`xi = 3`, `zeta = 2` et `4`). Le script définit donc la
**fissure principale** comme le plus court chemin, dans le graphe, de la pointe de l'entaille (option `--pointe`,
défaut (1.181 mm, 8 mm)) à l'extrémité de la même fissure la plus éloignée de cette pointe.

| Carte | Définition |
|---|---|
| tortuosité | longueur de la fissure principale / distance entre ses extrémités (1 = fissure droite) |
| part des ramifications | 1 − longueur de la fissure principale / longueur rompue totale (toutes fissures) |
| vitesse de la fissure principale | pente de sa longueur rompue en fonction du temps, ajustée entre 20 % et 60 % de sa longueur finale (m/s) ; au-delà de 60 %, certaines fissures s'arrêtent avant le bord |
| délai d'amorçage | de la première rupture à 20 % de la longueur de la fissure principale (ms) |

La grille `courbes_fissure_principale.png` donne la longueur de la fissure principale en fonction du temps depuis
la première rupture (mêmes axes pour tous les cas). `facies_crackwatch.png` dessine, pour chaque cas, le graphe
des faces rompues dans la configuration initiale (`conf0`) : fissure principale en bleu foncé, ramifications et
fissures secondaires en orange. Ce dessin remplace les `crack_N.svg` de `crackwatch`, avec la fissure principale
définie depuis l'entaille. Le cadrage vertical, commun à tous les cas, couvre la fissure issue de l'entaille et ses
ramifications : de petites fissures secondaires près des mors peuvent sortir du cadre. La vitesse ainsi mesurée varie de façon régulière avec `xi` et
`zeta`, contrairement à la vitesse moyenne de `plot_breaks.py`. Les délais d'amorçage longs (`xi = 2` et `3`,
`zeta ≥ 5`) sont réels : la fissure s'arrête après les premières ruptures en pointe d'entaille avant de repartir.

Les scripts sont ceux de la première étude, adaptés aux nouvelles sorties :

- la période de `top.txt` / `bottom.txt` est lue dans l'input de chaque cas (`define TCAPT`) ;
- les images `see2` montrent les champs et le faciès sans les parois des cellules (`show_cells=0`,
  `show_contours=0`), avec le cadre de l'échantillon (`show_bbox=1`, option ajoutée à `see2`) ;
- les conf rendus par `see2` sont choisis d'après leur temps : dernier conf (faciès de rupture), dernier conf avant
  le pic de force (champs), conf sauvegardé à la première rupture (directions principales « à l'initiation »).

Les grandeurs, les cartes et les options (`--vue`, `--echelle par_cas`…) sont décrites dans
`etude-param-model/README.md`, section « Dépouillement ».
