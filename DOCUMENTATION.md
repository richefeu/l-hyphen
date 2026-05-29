# Documentation d'utilisation de l-hyphen

## Dépendances

- **C++17** — compilateur compatible (g++, clang)
- **toofus** — bibliothèque utilitaire (clonée automatiquement)
- **GLFW3** — pour la visualisation 3D avec `see2` (optionnel)
- **OpenGL** — pour le rendu graphique
- **Typst** — pour recompiler le cheatsheet (optionnel, `brew install typst`)

## Compilation

```bash
make          # compile run et see2
make clean    # supprime les objets et les exécutables
make clean+   # supprime aussi le dossier toofus (cloné automatiquement)
```

La dépendance [toofus](https://github.com/richefeu/toofus) est clonée automatiquement si elle est absente.

---

## Exécutables

### `run` — lancer une simulation

```bash
./run <fichier_entrée>
```

Lit le fichier d'entrée, construit la configuration initiale et intègre les équations du mouvement.  
Les fichiers de sortie (`conf*`, `sample*.svg`) sont écrits dans le dossier courant.

### `see2` — visualiser une configuration

```bash
./see2                  # charge conf0 par défaut
./see2 <numéro>         # charge confN
./see2 <nom_de_fichier> # charge un fichier conf nommé
```

---

## Format du fichier d'entrée

Le fichier d'entrée est un fichier texte. Les lignes commençant par `#`, `/` ou `!` sont des commentaires. Les valeurs numériques acceptent des expressions arithmétiques et des constantes nommées (voir `define`).

### Constantes nommées

```
define <NOM> <valeur>
```
Définit une constante réutilisable dans les expressions suivantes.

---

### Paramètres de temps

| Commande | Description |
|----------|-------------|
| `dt <valeur>` | Pas de temps |
| `nstep <valeur>` | Nombre de pas à simuler |
| `t <valeur>` | Temps initial |
| `cyclicVelPeriod <valeur>` | Durée d'un cycle (inversion de la vitesse de chargement) |

---

### Dissipation

| Commande | Description |
|----------|-------------|
| `numericalDissipation <valeur>` | Coefficient de dissipation numérique |
| `globalViscosity <valeur>` | Viscosité globale |
| `setCellWallDampingRates <alpha_s> <alpha_b>` | Taux d'amortissement axial (`alpha_s`) et en flexion (`alpha_b`) |
| `setCellWallDampings <nu_s> <nu_b>` | Coefficients d'amortissement directs (axial et flexion) |

---

### Forces volumiques

```
gravity <gx> <gy>
```

---

### Détection de voisins

| Commande | Description |
|----------|-------------|
| `distVerlet <valeur>` | Distance de la zone tampon de Verlet |
| `nstepPeriodVerlet <valeur>` | Fréquence de mise à jour de la liste de voisins (en pas) |
| `linkCells <lx> <ly>` | Active l'algorithme link-cells (plus rapide pour grands systèmes) ; `lx`/`ly` sont les dimensions des cellules |

Sans `linkCells`, la recherche est en O(N²).

---

### Contact entre cellules

| Commande | Description |
|----------|-------------|
| `kn <valeur>` | Raideur normale de contact |
| `kt <valeur>` | Raideur tangentielle de contact |
| `viscn <valeur>` | Viscosité normale au contact |
| `mu <valeur>` | Coefficient de frottement |
| `fadh <valeur>` | Force d'adhésion normale |
| `adaptativeStiffness <0|1>` | Raideur adaptative (0 = désactivée, 1 = activée) |

---

### Cohésion (glue)

```
glue <dist>
```
Crée des points de colle entre nœuds distants de moins de `<dist>` (modèle à énergie critique).

```
GcGlue <dist>
distGcGlue <dist>
```
Idem avec le modèle basé sur l'énergie de rupture Gc.

```
setGlueSameProperties <kn_coh> <kt_coh> <fn_coh_max> <ft_coh_max> <yieldPower>
```
Définit les propriétés de cohésion uniformes pour toutes les zones collées (modèle force max).

```
setGcGlueSameProperties <kn_coh> <kt_coh> <Gc>
```
Idem pour le modèle Gc.

---

### Construction de l'échantillon

#### Lire un fichier de nœuds préprocessé

```
readNodeFile <fichier> <barWidth> <Kn> <Kr> <Mz_max> <p_int>
```

Le fichier contient des lignes `x y id_cellule`. Mettre `barWidth` négatif pour une estimation automatique (= moitié de la distance minimale entre nœuds de cellules différentes).

#### Ajouter un polygone régulier

```
addRegularPolygonalCell <nbFaces> <x> <y> <rot> <Rext> <barWidth> <kn> <kr> <mz_max>
```

#### Ajouter un arrangement mur de briques

```
addSquareBrickWallCells <nx> <ny> <hdist> <xleft> <ybottom> <barWidth> <kn> <kr> <mz_max>
```

#### Ajouter une ligne multi-barres

```
addMultiLine <xo> <yo> <xe> <ye> <barWidth> <n> <kn> <kr> <mz_max>
```

---

### Masses

| Commande | Description |
|----------|-------------|
| `setNodeMasses <masse>` | Masse identique pour tous les nœuds |
| `setCellMasses <masse>` | Masse totale par cellule (distribuée sur les nœuds) |
| `setCellWallDensities <rho> <épaisseur>` | Masse des nœuds à partir de la densité des parois |
| `setCellDensities <rho> <épaisseur>` | Idem en incluant la contribution surfacique |

---

### Contenu interne des cellules (pression)

| Commande | Description |
|----------|-------------|
| `cellContent <0\|1\|2>` | 0 = vide, 1 = PV constant (gaz), 2 = PV élastique (liquide) |
| `compressFactor <valeur>` | Facteur de compressibilité pour le mode élastique |
| `setCellInternalPressure <id_cellule> <p>` | Pression initiale d'une cellule |
| `setCellAsOpen <id_cellule>` | Marque la cellule comme ouverte (non fermée) |
| `setClose <id_cellule>` | Marque la cellule comme fermée |

---

### Conditions aux limites (contrôle des nœuds)

Le mode de contrôle : `0` = vitesse imposée, `1` = force imposée.

```
setNodeControl <id_cellule> <id_nœud> <xmode> <xvalue> <ymode> <yvalue>
```

```
setNodeControlInBox <xmin> <xmax> <ymin> <ymax> <xmode> <xvalue> <ymode> <yvalue>
```
Applique le contrôle à tous les nœuds dans la zone rectangulaire.

```
setCellControl <id_cellule> <xmode> <xvalue> <ymode> <yvalue>
```
Applique le contrôle à tous les nœuds de la cellule.

---

### Sorties

| Commande | Description |
|----------|-------------|
| `isvg <valeur>` | Index de départ des fichiers SVG |
| `nstepPeriodSVG <valeur>` | Fréquence de sauvegarde SVG (en pas ; 0 = désactivé) |
| `iconf <valeur>` | Index de départ des fichiers de configuration |
| `nstepPeriodConf <valeur>` | Fréquence de sauvegarde des configurations (en pas ; 0 = désactivé) |
| `nstepPeriodRecord <valeur>` | Fréquence d'enregistrement des données scalaires |
| `captureNodes <fichier> <xmin> <xmax> <ymin> <ymax>` | Enregistre les positions des nœuds dans la zone en post-traitement |
| `followCell <id_cellule>` | Suit une cellule particulière (données de suivi) |
| `findDisplayArea <facteur>` | Calcule automatiquement les limites d'affichage (facteur ≥ 1) |
| `limits <xmin> <xmax> <ymin> <ymax>` | Définit manuellement les limites d'affichage |

---

### Divers

| Commande | Description |
|----------|-------------|
| `nbThreads <n>` | Nombre de threads OpenMP (si compilé avec OpenMP) |
| `reorder <0\|1>` | Active/désactive le ré-ordonnancement des nœuds à la lecture |

---

## Visualiseur `see2` — raccourcis clavier

| Touche | Action |
|--------|--------|
| `a` | Afficher/masquer les zones de contrôle |
| `b` | Coloriser les barres selon l'effort axial |
| `c` | Afficher/masquer les cellules |
| `f` | Afficher/masquer les forces de contact |
| `g` | Afficher/masquer les points de colle |
| `h` | Afficher l'aide dans le terminal |
| `n` | Afficher/masquer les contours des cellules |
| `v` | Afficher/masquer les nœuds (points) |
| `p` | Afficher/masquer la pression interne |
| `q` | Quitter |
| `s` / `S` | Réduire / augmenter l'échelle des vecteurs force |
| `z` / `Z` | Zoom avant / arrière |
| `→` | Charger la configuration suivante (`confN+1`) |
| `←` | Charger la configuration précédente (`confN-1`) |
| `=` | Recadrer la vue sur l'ensemble de la scène |
| `↑` / `↓` | Augmenter / réduire le nombre de lignes de texte |

**Navigation à la souris :**
- Clic gauche + glisser → rotation/déplacement de la vue
- Maj + clic gauche + glisser → panoramique
- Clic milieu + glisser → zoom

---

### Modèles de contenu cellulaire

| Valeur | Modèle | Description |
|--------|--------|-------------|
| `0` | `CELL_EMPTY` | Cellule vide, sans modèle de pression |
| `1` | `CELL_ELASTIC_PV` | Compression élastique : $p = K \cdot \Omega / \Omega_0$ |
| `2` | `CELL_RIGID` | Cellule rigide (volume constant, non déformable) |

---

### Diagnostic et analyse

```bash
./run <fichier_entrée>
```

À chaque lancement, `run` affiche automatiquement un **diagnostic** contenant :
- **Géométrie** : nombre de cellules, nœuds, barres, voisins
- **Masses et raideurs** : gammes de valeurs
- **Stabilité** : rapport critique dt_crit/dt pour contact et cohésion
- **Verlet** : paramètres de la peau de détection
- **Sorties** : fréquences de sauvegarde

Le rapport complet est sauvegardé dans `diagnostic.txt`.

---

### Support du préprocesseur

Les fichiers d'entrée supportent :
- **Expressions arithmétiques** : `1.5 * 2 + 3.14` 
- **Constantes nommées** : `define pi 3.14159; ... dt <pi>`
- **Commentaires** : lignes commençant par `#`, `/`, ou `!`

---

## Cheatsheet complet

Pour une référence rapide de toutes les commandes, voir [`cheatsheets/lhyphen_cheatsheet.pdf`](./cheatsheets/lhyphen_cheatsheet.pdf).

Compiler le cheatsheet :
```bash
cd cheatsheets && make build
```

---

## Fichiers générés

| Fichier | Description |
|---------|-------------|
| `conf0`, `conf1`, … | Configurations binaires lisibles par `see2` et `run` |
| `sample0000.svg`, … | Captures SVG de l'état du système |
| `diagnostic.txt` | Rapport de diagnostic (géométrie, stabilité, paramètres Verlet) |
| `breakHistory.txt` | Historique des ruptures de colle |
| `breakEvol.txt` | Évolution cumulée des ruptures |
| `lhyphen0000.txt` | Données scalaires enregistrées (fréquence `nstepPeriodRecord`) |

---

## Concepts clés

### Cellules et nœuds

- Une **cellule** est un polygone fermé (ou ouvert).
- Chaque cellule est composée de **nœuds** (sommets) reliés par des **barres** (arêtes).
- Les **nœuds** sont les points de calcul mécanique ; ils ont masse, vitesse, forces.
- Les **barres** connectent les nœuds ; elles transmettent efforts axiaux et moments de flexion.

### Détection de voisins (Verlet)

L'algorithme Verlet optimise la détection de contacts en maintenant une **liste de voisins** mise à jour à intervalles réguliers (`nstepPeriodVerlet`).

- **`distVerlet`** : épaisseur de la "peau" autour de chaque cellule.
- Si la distance entre deux cellules < (rayon1 + rayon2 + `distVerlet`), elles sont voisines.
- Entre deux mises à jour, les contacts supposés rester valides selon le déplacement maximal.

Pour grandes simulations (N > 500 nœuds), activez `linkCells` pour passer d'O(N²) à O(N).

### Cohésion (glue)

Deux modèles de rupture sont disponibles :

1. **Modèle force maximale** (`glue`, `setGlueSameProperties`) : rupture si l'effort dépasse un seuil.
2. **Modèle énergie critique** (`GcGlue`, `setGcGlueSameProperties`) : rupture si l'énergie élastique dépasse Gc (énergie par unité de longueur).

---

## Conseils d'utilisation

### Choix du pas de temps

Le diagnostic indique le ratio `dt_crit / dt` pour contact et cohésion :
- **≥ 50** : dt peut être augmenté
- **≥ 20** : stable et sûr
- **10–20** : proche de la limite, vérifier la stabilité
- **< 10** : instabilité imminente

### Performance

| Amélioration | Action |
|-------------|--------|
| Réduire le temps CPU | Augmenter `nstepPeriodVerlet` ou activer `linkCells` |
| Réduire la mémoire | Diminuer les fréquences de sauvegarde (`nstepPeriodSVG`, `nstepPeriodConf`) |
| Améliorer la précision | Réduire `dt` (mais réduit la stabilité) |
| Accélérer l'intégration | Augmenter `numericalDissipation` (mais perd la physique) |

### Troubleshooting

| Problème | Cause probable | Solution |
|----------|---|---|
| Instabilité explosif | dt trop grand | Réduire dt, vérifier `dt_crit/dt` dans le diagnostic |
| Contacts manqués | `distVerlet` trop petit ou `nstepPeriodVerlet` trop grand | Augmenter `distVerlet` ou réduire `nstepPeriodVerlet` |
| Cellules pénétrantes | `kn` trop faible ou `adaptativeStiffness = 0` | Augmenter `kn` ou activer `adaptativeStiffness = 1` |
| Sorties bloquées | Disque plein ou permissions insuffisantes | Vérifier l'espace disque ou augmenter les fréquences |

---

## Exemple minimal de fichier d'entrée

```
# Paramètres de simulation
dt          1e-6
nstep       100000

numericalDissipation  0.0
globalViscosity       2.0

gravity     0  -9.81

distVerlet          0.05
nstepPeriodVerlet   50

# Sorties
iconf             0
nstepPeriodConf   10000
isvg              0
nstepPeriodSVG    10000

# Contact
kn   1e4
kt   1e4
mu   0.3

# Échantillon : 6 cellules hexagonales sur grille triangulaire
#   addRegularPolygonalCell nbFaces x y rot Rext barWidth kn kr mz_max
addRegularPolygonalCell  6  0.0  0.0  0.0  0.5  0.05  1e4  500  1e9

setNodeMasses  0.1

# Conditions aux limites (vitesse nulle en bas)
setNodeControlInBox  -1.0  1.0  -0.1  0.01  0  0.0  0  0.0

findDisplayArea  1.1
```

---

## Ressources supplémentaires

- **GitHub** : https://github.com/richefeu/l-hyphen
- **Cheatsheet complet** : `cheatsheets/lhyphen_cheatsheet.pdf`
- **Bibliothèque toofus** : https://github.com/richefeu/toofus
- **Exemples** : dossier `examples/`

---

## Glossaire rapide

| Terme | Sens |
|-------|------|
| **Barre** | Élément structurel connectant deux nœuds ; transmet axial et flexion |
| **Cellule** | Polygone fermé (ou ouvert) composé de nœuds et barres |
| **Cohésion** | Liaison glued entre deux cellules, avec critère de rupture |
| **dt_crit** | Pas de temps critique pour stabilité (Courant-Friedrichs-Lewy) |
| **Glue** | Point de collage entre cellules voisines |
| **kn, kt** | Raideurs : normales et tangentielles pour contact |
| **kr** | Raideur de flexion (moment) au nœud |
| **Nœud** | Point de calcul mécanique (masse, vitesse, accélération) |
| **Rayon (radius)** | Demi-épaisseur fictive des barres pour calculs géométriques |
| **Verlet** | Algorithme de détection de voisins avec peau tampon |

---

## Historique des modifications

- **2026-05-29** : Amélioration du diagnostic (suppression dt_crit barres/flexion, meilleure lisibilité, sauvegarde en fichier)
- **2026-05-29** : Complétion du cheatsheet et création du Makefile pour sa compilation
- **2026-05-29** : Amélioration et complétion de cette documentation
