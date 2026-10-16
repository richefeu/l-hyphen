# nodeGen — génération et visualisation de nodeFiles

Deux petits outils C++17 autonomes (aucune dépendance) :

- **`nodegen`** génère un nodeFile (`x y idCellule`, lu par `readNodeFile`) : pavage de Voronoï d'un
  rectangle, avec une pré-fissure rectiligne nette et sans barres trop courtes ;
- **`nodeview`** produit un SVG d'analyse de n'importe quel nodeFile (parois collées ou libres, barres courtes).

```bash
cd nodeGen && make
cd example && ../nodegen params.txt          # -> nodefile.txt + nodefile.svg
../nodeview nodefile.txt                     # -> nodefile.svg
../nodeview nodefile.txt -box 0.0040 0.0047 0.0097 0.0103 -nodes -o tip.svg   # zoom pointe de fissure
```

Le SVG s'ouvre dans un navigateur (zoom illimité, c'est un format vectoriel).

Un tutoriel détaillé (utilisation, algorithmes, validation, améliorations possibles) est dans
`doc/tutorial_nodegen.pdf` ; ses figures se régénèrent avec `doc/make_figures.sh`.

## nodegen

### Principe

1. **Germes** (`lattice`) :
   - `hex` (défaut) : réseau hexagonal de même densité, dont chaque germe est déplacé aléatoirement d'au plus
     `disorder` × le pas du réseau. Le pavage est alors dominé par des hexagones, avec des pentagones et quelques
     heptagones d'autant plus nombreux que `disorder` est grand. En présence d'une pré-fissure, le réseau est
     orienté selon celle-ci (rangées parallèles et symétriques par rapport à la fissure).
   - `random` : `N = Lx·Ly / (π·cellSize²/4)` points tirés avec une distance minimale (0.7 × espacement moyen) ;
     pavage plus désordonné (environ 20 % de pentagones, 60 % d'hexagones et 20 % d'heptagones).
2. **Régularisation de Lloyd** (`lloydIterations`) : chaque germe est déplacé au centroïde de sa cellule.
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

### Paramètres (`params.txt`)

| Mot-clé | Défaut | Description |
|---|---|---|
| `xmin`, `ymin` | 0, 0 | coin inférieur gauche du domaine |
| `Lx`, `Ly` | 0.01 | dimensions du domaine |
| `cellSize` | 3.5e-4 | diamètre équivalent visé des cellules |
| `barWidth` | 2e-6 | largeur des barres (= `barWidth` de `readNodeFile`) |
| `lattice` | hex | `hex` : réseau hexagonal perturbé ; `random` : tirage aléatoire |
| `disorder` | 0.3 | perturbation des germes du réseau `hex`, en fraction du pas (0 : hexagones parfaits) |
| `lloydIterations` | 10 | itérations de régularisation de Lloyd |
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
| `grips h v` | — | forme historique de `gripHeight h` + `pullVelocity v` |
| `nodeMass m` | — | masse des nœuds, pour estimer le pas de temps critique de flexion |
| `input` | nodegen_input.txt | lignes d'input l-hyphen produites |
| `inputTemplate modèle [sortie]` | —, input.txt | input complet produit à partir d'un input existant |
| `output` | nodefile.txt | nodeFile produit |
| `svg` | nodefile.svg | vue d'analyse produite |

Les lignes commençant par `#` (ou la fin de ligne après `#`) sont des commentaires.

## nodeview

```
nodeview nodefile.txt [-w barWidth] [-g distGlue] [-r ratio] [-box x0 x1 y0 y1] [-nodes] [-bare] [-px largeur] [-o out.svg]
```

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
