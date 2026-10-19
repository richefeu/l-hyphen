# nodeGen — génération et visualisation de nodeFiles

Deux petits outils C++17 autonomes (aucune dépendance) :

- **`nodegen`** génère un nodeFile (`x y idCellule`, lu par `readNodeFile`) : pavage de Voronoï d'un
  rectangle, avec une pré-fissure rectiligne nette et sans barres trop courtes, et éventuellement un tissu
  végétal organisé autour de faisceaux vasculaires (bundles, voir [Tissu végétal](#tissu-végétal--bundles-et-attraction-de-mahalanobis)) ;
- **`nodeview`** produit un SVG d'analyse de n'importe quel nodeFile (parois collées ou libres, barres courtes).

```bash
make                                         # dans nodeGen/
cd exampleTissu && ../nodegen params.txt     # tissu avec bundles -> nodefile.txt, nodefile.svg, bundles.txt
cd ../example && ../nodegen params.txt       # échantillon pré-fissuré -> nodefile.txt + nodefile.svg
../nodeview nodefile.txt                     # -> nodefile.svg
../nodeview nodefile.txt -box 0.0040 0.0047 0.0097 0.0103 -nodes -o tip.svg   # zoom pointe de fissure
```

Le SVG s'ouvre dans un navigateur (zoom illimité, c'est un format vectoriel).

Un tutoriel détaillé (utilisation, tissus végétaux, algorithmes, validation, limites) est dans
`doc/tutorial_nodegen.pdf` ; ses figures et balayages se régénèrent avec `doc/make_figures.sh`.

## nodegen

### Principe

1. **Germes** (`lattice`) :
   - `hex` (défaut) : réseau hexagonal de même densité, dont chaque germe est déplacé aléatoirement d'au plus
     `disorder` × le pas du réseau. Le pavage est alors dominé par des hexagones, avec des pentagones et quelques
     heptagones d'autant plus nombreux que `disorder` est grand. En présence d'une pré-fissure, le réseau est
     orienté selon celle-ci (rangées parallèles et symétriques par rapport à la fissure).
   - `square` : grille carrée de même densité (même perturbation `disorder`, même orientation selon la fissure) ;
     avec `disorder 0`, c'est la grille carrée exacte utilisée dans l'article pour la comparaison aux coupes
     tomographiques.
   - `random` (ou `poisson`) : `N = Lx·Ly / (π·cellSize²/4)` points tirés avec une distance minimale
     (0.7 × espacement moyen) ; pavage plus désordonné (environ 20 % de pentagones, 60 % d'hexagones et 20 % d'heptagones).
   - `uniform` : `N` points tirés uniformément, sans distance minimale (processus de Poisson).
2. **Régularisation de Lloyd** (`lloydIterations`) : chaque germe est déplacé au centroïde de sa cellule ; avec
   des bundles et `alpha > 0`, ce centroïde est en plus attiré vers le bundle le plus proche (voir plus bas).
3. **Pavage de Voronoï** découpé au rectangle (découpe par demi-plans, cellules convexes).
4. **Arêtes courtes et jonctions à 3 branches.** Une jonction où se rejoignent plus de 3 parois est
   mécaniquement indésirable : tous les sommets du pavage ont 3 branches (bords compris), sauf les extrémités de
   la pré-fissure. Les arêtes courtes ne sont donc **jamais fusionnées**, puisque la fusion crée un sommet à
   4 branches. Une arête est courte si elle mesure moins de `minEdgeFactor` × `lRef`, avec
   `lRef = sqrt(Lx·Ly / Ncell)` la longueur caractéristique des cellules : c'est le côté du carré de même aire
   qu'une cellule (`Ncell` est le nombre de cellules du pavage). `lRef` ne dépend que de la densité de cellules,
   alors que la longueur moyenne des arêtes `l_moy` (environ 0.62 `lRef` pour des hexagones) dépend aussi de leur
   forme, donc du traitement lui-même. Les arêtes courtes sont traitées :
   - par une **transition T1** si les deux cellules qui la bordent gardent au moins `minSides` côtés : l'arête
     est remplacée par une arête transversale de longueur `t1Margin` × seuil, et ses deux cellules cessent
     d'être voisines tandis que les deux autres le deviennent ;
   - sinon par **étirement** : ses deux sommets sont écartés le long de l'arête jusqu'à un peu plus que le
     seuil, ce qui conserve le nombre de côtés. Les sommets sur les bords du domaine ou sur la fissure restent
     sur leur ligne, et un écartement qui rendrait une cellule non convexe est annulé ;
   - en dernier recours, un **amas** d'arêtes courtes voisines est contracté en un sommet, aussitôt éclaté en
     sommets à 3 branches, sans faire passer une cellule sous `minSides` côtés.

   Les sommets à 4 branches ou plus du pavage de Voronoï (sommets quasi confondus, ligne de symétrie prolongeant
   la fissure, voir plus bas) sont éclatés en sommets à 3 branches. Les barres très courtes, qui imposent un pas
   de temps minuscule en flexion (`kr/l²`), sont ainsi évitées à la source (voir
   `Frank-ouverture-fissure/note_barres_courtes`).

   `nodegen` vérifie le résultat. Il affiche le nombre de sommets par nombre de branches et **s'arrête avec un
   code d'erreur (2)** s'il reste un sommet à 4 branches hors des extrémités de la fissure, ou une arête sous le
   seuil. Les fichiers sont tout de même écrits pour pouvoir être examinés. Sur la géométrie de l'exemple, 20
   graines par réglage :

   | `minEdgeFactor` | `hex`, `disorder 0.3` | `hex`, `disorder 0.6` | `random`, `lloydIterations 10` |
   |---|---|---|---|
   | 0.22 | 0 échec | 0 échec | 0 échec |
   | 0.25 | — | 1 échec | 0 échec |
   | 0.28 | — | 2 échecs | 0 échec |
   | 0.33 | 0 échec | 3 échecs | 3 échecs |

   Changer de graine suffit en général.
5. **Décalage des parois** de `barWidth/2` vers l'intérieur (même construction que cellPrepro) : deux cellules
   voisines sont exactement à `barWidth` l'une de l'autre, condition du collage par `glue` / `distGcGlue`.

Répartition obtenue sur la géométrie de l'exemple (`minEdgeFactor 0.22`, `minSides 5`, cellules intérieures) :

| Réglage | 4 côtés | 5 côtés | 6 côtés | 7 côtés | 8 côtés |
|---|---|---|---|---|---|
| `lattice random`, `lloydIterations 10` | 0 % | 22 % | 57 % | 20 % | 0.7 % |
| `lattice hex`, `disorder 0.2`, `lloydIterations 5` | 0.1 % | 2.3 % | 96 % | 1.4 % | 0.1 % |
| `lattice hex`, `disorder 0.3`, `lloydIterations 5` | 0 % | 2.9 % | 95 % | 1.7 % | 0.1 % |
| `lattice hex`, `disorder 0.4`, `lloydIterations 5` | 0.1 % | 7.6 % | 85 % | 6.7 % | 0 % |
| `lattice hex`, `disorder 0.6`, `lloydIterations 5` | 0.1 % | 26 % | 51 % | 21 % | 1.9 % |

Avec un seuil plus élevé (`minEdgeFactor 0.33`, réglage de `example/params.txt`), les transitions T1 sont plus
nombreuses : avec `disorder 0.3`, 92 % d'hexagones ; avec `disorder 0.6`, 48 %.

Une transition T1 retire un côté à deux cellules et en ajoute un à deux autres. Par rapport à l'ancienne fusion,
qui ne faisait que retirer des côtés, il y a donc plus d'heptagones et moins de pentagones.

Les cellules de bord, coupées par le rectangle, ont naturellement moins de côtés ; `nodegen` affiche la
répartition sur toutes les cellules et sur les cellules intérieures.

### Pré-fissure

La pré-fissure est un segment `crack x0 y0 x1 y1`. Les germes situés dans une bande autour de ce segment sont
placés par paires symétriques par rapport à la ligne de fissure ; par symétrie, la ligne est alors exactement
formée d'arêtes du pavage. Aucune cellule n'est coupée ni retirée : la fissure est rectiligne et les cellules
de part et d'autre sont complètes.

Le long de la fissure, les parois sont décalées de `barWidth/2 + opening/2` : l'écart entre les deux lèvres
vaut `barWidth + opening`, supérieur à `barWidth + distGcGlue`, donc ces parois ne sont pas collées.

**Pointe de fissure.** La symétrie est prolongée au-delà de la pointe : la ligne de fissure s'y continue par
une interface **collée** entre deux cellules. La fissure se termine donc sur un sommet où se rejoignent quatre
cellules, en face d'une interface, et non contre une cellule. C'est la seule jonction à 4 branches conservée.
La pointe est le sommet de la ligne le plus proche de `(x1, y1)` ; sa position effective est affichée.

Au-delà de la pointe, chaque sommet de la ligne de symétrie a 4 branches (deux cellules de chaque côté). Il est
éclaté en deux sommets à 3 branches qui rendent voisines deux cellules en diagonale : l'interface collée qui
prolonge la fissure n'est donc pas rectiligne mais en **zigzag**. Le long de la fissure ouverte, les sommets
gardent leur structure : les lèvres sont libres, et il n'y a mécaniquement que 3 branches de chaque côté.

Une ligne droite d'arêtes n'existe pas dans un pavage hexagonal : la rangée de pentagones qui borde la fissure
doit se raccorder au réseau régulier, ce qui crée un défaut (dislocation : quelques pentagones et heptagones)
au bout de l'interface prolongeant la fissure. `tipInterface` fixe la longueur de cette interface, en nombre de
cellules, donc la distance entre la pointe et ce défaut. Au-delà, le réseau est régulier : la fissure n'est pas
prolongée par une interface rectiligne qui serait un chemin de rupture privilégié.

La fissure peut partir d'un bord (fissure débouchante) ou être entièrement intérieure, avec une orientation
quelconque.

### Tissu végétal : bundles et attraction de Mahalanobis

Cette partie reprend le modèle du draft `TissuGenerationDraft/` (Trivino *et al.*, *A continuous description of
plant cellular tissues from variational Voronoi morphogenesis*). Un tissu végétal y est vu comme deux populations
de cellules : les cellules ordinaires, qui remplissent l'espace, et les **bundles** (faisceaux vasculaires), qui
apportent eau et nutriments. Le tissu est le minimiseur de la fonctionnelle

  E(P, B) = E_CVL(P) − α E_att(P, B),

où E_CVL = Σᵢ ∫_{Vᵢ} ‖x − pᵢ‖² dx est l'énergie de Lloyd (compacité, tendance à l'hexagone régulier) et
E_att = Σᵢ ‖pᵢ − b_{j(i)}‖²_M attire chaque cellule vers son bundle le plus proche au sens de la distance de
Mahalanobis d(c, b) = √((c − b)ᵀ Σ⁻¹ (c − b)).

**Algorithme.** Une fraction `bundleFraction` des germes est tirée au hasard comme bundles (après la construction
de l'état initial, avec la graine `seed`), plus le germe le plus proche de chaque point `bundle x y`. Ensuite, à
chacune des `lloydIterations` itérations, pour chaque germe ordinaire :

1. calcul du centroïde c de sa cellule de Voronoï (terme E_CVL) ;
2. recherche du bundle b le plus proche de c au sens de Mahalanobis ;
3. nouveau germe p' = (c + α b) / (1 + α) : c est déplacé vers b d'une fraction α/(1+α) du segment [c, b].

Tous les germes sont mis à jour simultanément, avec les positions des bundles du début de l'itération. Les bundles
eux-mêmes suivent `bundleMotion` : `lloyd` (défaut, ils vont au centroïde de leur cellule ; B est optimisé avec P),
`fixed` (immobiles) ou `attract` (attirés comme les autres cellules, vers le bundle le plus proche autre qu'eux).
Avec `alpha 0` (défaut), on retrouve exactement la régularisation de Lloyd ; sans bundle, `alpha` est sans effet.

**Les deux paramètres de forme.**

- `alpha` règle surtout l'**hétérogénéité des aires** : les cellules proches d'un bundle rétrécissent (l'annexe
  du draft montre que l'aire décroît strictement quand le germe est attiré). Ordres de grandeur de l'article :
  α entre 0 et 0.2, référence α = 0.05.
- l'anisotropie `K` de Σ règle l'**organisation spatiale** : `K = 1` (Σ = I) donne une distance euclidienne et des
  amas de petites cellules autour de chaque bundle ; `K ≫ 1` (Σ = diag(1, K)) rend la distance presque
  insensible aux écarts selon l'axe des colonnes : un bundle attire toute la colonne de cellules de même abscisse
  et le tissu s'organise en bandes (K = 1000 dans l'article, K = 10000 pour les coupes tomographiques).
  `anisotropy K angle` oriente les colonnes (angle en degrés depuis l'axe x, 90 par défaut, soit Σ = diag(1, K)) ;
  `sigma sxx sxy syy` donne une matrice quelconque (définie positive). Seule la forme de Σ compte, pas son échelle.

Recettes de l'article (domaine et `cellSize` à adapter) :

| Tissu | Réglage |
|---|---|
| hexagonal pur | `alpha 0` (Lloyd seul) |
| référence (amas autour des bundles) | `lattice random`, `bundleFraction 0.05`, `alpha 0.05`, `anisotropy 1` |
| colonnes | idem avec `anisotropy 1000` |
| coupe tomographique | `lattice square`, `disorder 0`, `lloydIterations 5`, `bundleFraction 0.05`, `alpha 0.02`, `anisotropy 10000` |

Une forte anisotropie crée des cellules très inégales et davantage d'arêtes courtes : vérifier la ligne
« arêtes étirées … encore trop courtes » et `l_min/l_moy` dans la sortie (pas de temps critique de flexion).

**Sorties.**

- `bundleOutput` (défaut `bundles.txt`) : identifiants des cellules bundles dans le nodeFile. `readNodeFile`
  numérote les cellules dans l'ordre des identifiants : ce sont aussi les indices des cellules dans l-hyphen si
  le nodeFile est la première source de cellules de l'input (utile pour leur donner des propriétés propres). Le
  fichier n'est écrit que s'il y a des bundles ;
- le SVG colorie les bundles (`nodeview nodefile.txt -fill bundles.txt` fait de même pour un nodeFile existant) ;
- les **métriques** de l'article sont affichées : indice de forme moyen ⟨SI⟩ = ⟨4πA/P²⟩ (compacité, 0.907 pour
  l'hexagone régulier), coefficient de variation des aires σ_A/⟨A⟩ (hétérogénéité), distance moyenne au plus
  proche voisin ⟨d_NN⟩ (entre centroïdes) et distance de Mahalanobis moyenne au bundle le plus proche ⟨d_M⟩,
  ces deux dernières rapportées à l'espacement visé √(π cellSize²/4) ; aire, SI et d_NN sont aussi donnés
  séparément pour les cellules ordinaires et les bundles. Elles sont calculées sur le pavage avant le décalage
  des parois, sur toutes les cellules (bords compris) ;
- `morphology fichier` ajoute une ligne de ces métriques au fichier (en-tête écrit à la création) : un balayage
  en α, K ou nombre d'itérations se fait en relançant `nodegen` avec des paramètres différents et le même
  fichier `morphology`.

**Nombre d'itérations.** L'itération ne converge pas vers un état stationnaire : à α et K fixés, les cellules
continuent de se resserrer autour des bundles. À K = 1, σ_A/⟨A⟩ croît avec `lloydIterations` sans saturer
(α = 0.05 : 0.32 après 5 itérations, 0.64 après 50, 0.83 après 200) ; à K = 1000, elle passe par un maximum vers
30–50 itérations puis décroît. `lloydIterations` est donc un troisième paramètre de
forme ; deux tissus ne se comparent qu'à nombre d'itérations égal.

Les bundles sont compatibles avec la pré-fissure : la symétrie par rapport à la fissure est rétablie à chaque
itération et un bundle situé dans la bande symétrique devient une paire de bundles symétriques.

### Paramètres (`params.txt`)

| Mot-clé | Défaut | Description |
|---|---|---|
| `xmin`, `ymin` | 0, 0 | coin inférieur gauche du domaine |
| `Lx`, `Ly` | 0.01 | dimensions du domaine |
| `cellSize` | 3.5e-4 | diamètre équivalent visé des cellules |
| `barWidth` | 2e-6 | largeur des barres (= `barWidth` de `readNodeFile`) |
| `lattice` | hex | état initial : `hex` réseau hexagonal perturbé ; `square` grille carrée perturbée ; `random` (ou `poisson`) tirage avec distance minimale ; `uniform` tirage uniforme |
| `disorder` | 0.3 | perturbation des germes des réseaux `hex` et `square`, en fraction du pas (0 : réseau parfait) |
| `lloydIterations` | 10 | itérations de régularisation de Lloyd (attraction vers les bundles comprise) |
| `bundleFraction f` | 0 | fraction des germes tirés au hasard comme bundles |
| `bundle x y` | — | bundle imposé : germe le plus proche du point (répétable) |
| `alpha` | 0 | intensité de l'attraction vers le bundle le plus proche, p' = (c + α b)/(1 + α) |
| `anisotropy K [angle]` | 1 90 | Σ d'anisotropie K, colonnes orientées à `angle` degrés de l'axe x (90 : Σ = diag(1, K)) |
| `sigma sxx sxy syy` | 1 0 1 | matrice de covariance Σ complète (distance de Mahalanobis) |
| `bundleMotion` | lloyd | déplacement des bundles : `lloyd`, `fixed` ou `attract` |
| `bundleOutput` | bundles.txt | identifiants des cellules bundles |
| `morphology fichier` | — | ajoute une ligne de métriques morphologiques au fichier |
| `minEdgeFactor` | 0.2 | seuil des arêtes courtes, en fraction de `lRef = sqrt(Lx·Ly / Ncell)` (0 : aucun traitement, déconseillé). L'ancien mot-clé `minEdgeRatio` (fraction de la longueur moyenne des arêtes) est refusé avec la valeur équivalente, environ 0.62 × `minEdgeRatio` |
| `t1Margin` | 1.2 | longueur de l'arête créée par une transition T1 ou un éclatement, en multiple du seuil |
| `minSides` | 5 | nombre de côtés minimal d'une cellule après une transition T1 (sinon : étirement) |
| `seed` | 1 | graine du générateur aléatoire |
| `crack x0 y0 x1 y1` | — | pré-fissure (optionnelle) |
| `tipInterface` | 1 | longueur de l'interface collée qui prolonge la fissure au-delà de la pointe (en cellules) |
| `opening` | 2e-6 | ouverture supplémentaire de la fissure (doit dépasser `distGlue`) |
| `distGlue` | 2e-7 | valeur de `distGcGlue` / `glue` de l'input (contrôle et SVG) |
| `cellProperties Kn Kr Mz_max p_int` | — | propriétés reportées dans la ligne `readNodeFile` générée |
| `glueProperties kn_coh kt_coh Gc` | — | propriétés reportées dans la ligne `setGcGlueSameProperties` générée |
| `gripRows n` | 1 | nombre de rangées de nœuds attrapées par chacun des mors bas et haut |
| `gripHeight h` | — | hauteur explicite des mors (utilisée si `gripRows` n'est pas donné) |
| `pullVelocity v` | 3e-4 | vitesse imposée aux mors (bas −v, haut +v) |
| `nodeMass m` | — | masse des nœuds, pour estimer le pas de temps critique de flexion |
| `input` | nodegen_input.txt | lignes d'input l-hyphen produites |
| `inputTemplate modèle [sortie]` | —, input.txt | input complet produit à partir d'un input existant |
| `output` | nodefile.txt | nodeFile produit |
| `svg` | nodefile.svg | vue d'analyse produite |

Les lignes commençant par `#` (ou la fin de ligne après `#`) sont des commentaires.

## nodeview

```
nodeview nodefile.txt [-w barWidth] [-g distGlue] [-r ratio] [-box x0 x1 y0 y1] [-nodes] [-bare] [-px largeur] [-o out.svg] [-fill ids.txt]
```

`-fill ids.txt` colorie les cellules dont l'identifiant est listé (par exemple `bundles.txt` de `nodegen`).

- **gris** : paroi collée (son milieu est à moins de `barWidth + distGlue` d'une barre d'une autre cellule,
  critère de création des liens cohésifs) ;
- **orange** : paroi libre (bord du domaine, pré-fissure, trou) ;
- **bleu** : barre plus courte que `ratio` × longueur moyenne. Dans le SVG produit par `nodegen`, ce sont les
  barres dont l'arête du pavage est sous le seuil `minEdgeFactor` × `lRef`, c'est-à-dire exactement les arêtes
  signalées en `ERREUR` par le contrôle final (deux barres par arête, une par cellule). Le décalage des parois
  raccourcit chaque barre d'environ `barWidth` : une barre à peine sous le seuil, dont l'arête est au seuil, n'est
  pas signalée ;
- **cercle rouge** : jonction à 4 branches ou plus, c'est-à-dire au moins 4 cellules en un sommet, ou 3 cellules
  sur un bord libre. Avec `nodegen`, seules les extrémités de la pré-fissure doivent apparaître. Les jonctions
  sont reconstruites à partir des parois collées : les lèvres libres d'une fissure ne sont pas rattachées entre
  elles.

`nodeview` affiche aussi le nombre de jonctions par nombre de cellules.

Sans `-w`, `barWidth` est estimé comme la distance minimale entre nœuds de cellules différentes. Le nombre de
parois collées affiché correspond au nombre de « Liens collés » de `diagnostic.txt`.

## Utilisation dans l-hyphen

`nodegen` écrit dans `input` (défaut `nodegen_input.txt`) les lignes d'input qui dépendent du maillage, à copier
dans le fichier d'input aux endroits indiqués :

- `readNodeFile` avec le bon `barWidth` (et `Kn Kr Mz_max p_int` si `cellProperties` est donné) ;
- les mors : deux `setNodeControlInBox` (bas et haut, libres en x, vitesse imposée en y) et les deux
  `captureNodes` correspondants (`bottom.txt`, `top.txt`), avec le nombre de nœuds capturés ;
- `distGcGlue` (et `setGcGlueSameProperties` si `glueProperties` est donné) ;
- en commentaire : statistiques du maillage, position effective de la pointe et, si `nodeMass` et `Kr` sont
  connus, le pas de temps critique de flexion `l_min·sqrt(m/Kr)`.

Une ligne dont une valeur n'est pas connue est écrite en commentaire avec des repères `<...>`. Exemple :

```
readNodeFile  nodefile.txt  2e-06  2000  0.000792863858442045  79286385.8442045  0
setNodeControlInBox  0.0032  0.0168  0.002  0.0024  1 0.0  0 -0.0003    # bas  : 312 noeuds
setNodeControlInBox  0.0032  0.0168  0.0176  0.018  1 0.0  0 0.0003    # haut : 316 noeuds
captureNodes  bottom.txt  0.0032  0.0168  0.002  0.0024
captureNodes  top.txt     0.0032  0.0168  0.0176  0.018
distGcGlue 2e-07
setGcGlueSameProperties  26373.6263736264  26373.6263736264  1.25751147880661e-06
```

### Mors : rangées de nœuds attrapées

Par défaut (`gripRows 1`), chaque mors attrape **une seule rangée de nœuds** : les nœuds du bord du domaine,
tous à `barWidth/2` du bord. Les rangées sont définies topologiquement : la rangée k est formée des sommets
situés à k−1 barres du bord. Une boîte de `setNodeControlInBox` ne sélectionne que par la hauteur ; `nodegen` la
fait s'arrêter à mi-distance entre le nœud le plus haut des `gripRows` premières rangées et le nœud suivant.

- `gripRows 1` : la rangée du bord est attrapée exactement (le nœud suivant est à environ un quart de cellule) ;
- `gripRows n` (n ≥ 2) : toutes les rangées demandées sont attrapées, mais un pavage de Voronoï n'a pas de rangées
  horizontales au-delà de la première ; la boîte attrape donc aussi quelques nœuds plus bas des rangées suivantes.
  Leur nombre est affiché (par exemple `253 noeuds (2 rangées + 45 noeuds des rangées suivantes)`) ;
- `gripHeight h` impose une hauteur de boîte.

### Input complet à partir d'un modèle

Avec `inputTemplate modele.txt input.txt`, `nodegen` recopie un input l-hyphen existant (par exemple celui d'une
étude précédente) dans `input.txt`, directement exécutable par `run`, en remplaçant les lignes qui dépendent du
maillage ; elles sont marquées `# [nodegen]` :

| Ligne du modèle | Traitement |
|---|---|
| `readNodeFile` | nouveau nodeFile (chemin relatif au répertoire de `input.txt`) et `barWidth` ; `Kn Kr Mz_max p_int` de `cellProperties`, sinon recopiés du modèle |
| `cleanShortBars` | commentée (inutile) |
| 2 × `setNodeControlInBox` | boîtes recalculées sur le nouveau domaine, de hauteur `gripRows` (ou `gripHeight`) ; modes et valeurs repris du modèle, sauf si `pullVelocity` est donné |
| 2 × `captureNodes` | boîtes recalculées, noms de fichiers conservés |
| `setGcGlueSameProperties` | remplacée si `glueProperties` est donné, sinon recopiée |
| `distGcGlue` / `GcGlue` / `glue` | recopiée ; `nodegen` vérifie qu'elle reste inférieure à `opening` |
| tout le reste | recopié tel quel |

Les mors bas et haut sont reconnus comme les deux boîtes de plus petit et de plus grand `ymin` ; si le modèle ne
contient pas exactement deux `setNodeControlInBox` (ou deux `captureNodes`), ces lignes sont recopiées sans
modification et un avertissement est affiché. `nodegen` affiche aussi le nombre de nœuds de chaque mors et, si la
masse des nœuds et `Kr` sont connus, le pas de temps critique de flexion. Il refuse d'écraser le modèle.

`cleanShortBars` n'est pas nécessaire avec un nodeFile produit par `nodegen`.
