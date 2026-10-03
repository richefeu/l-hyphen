#!/bin/bash
# Régénère les figures du tutoriel nodegen (nécessite inkscape et ImageMagick).
# Usage : cd nodeGen/doc && ./make_figures.sh
set -e
cd "$(dirname "$0")"
NG=../nodegen
NV=../nodeview
mkdir -p figures work
( cd .. && make -s )

png() { inkscape "$1" --export-type=png --export-filename="$2" -w "$3" > /dev/null 2>&1; }

# Petit domaine de démonstration (3 mm x 3 mm, pré-fissure de 1 mm depuis le bord gauche)
demo() { # nom puis paires clé=valeur
  local name=$1; shift
  cat > work/$name.params <<EOF
xmin 0
ymin 0
Lx 0.003
Ly 0.003
cellSize 3.53e-4
barWidth 2e-6
seed 3
crack 0 0.0015 0.001 0.0015
opening 2e-6
distGlue 2e-7
output work/$name.txt
svg work/$name.svg
input work/${name}_input.txt
EOF
  for kv in "$@"; do echo "${kv%%=*} ${kv#*=}" >> work/$name.params; done
  $NG work/$name.params > work/$name.log
  $NV work/$name.txt -bare -px 700 -o work/${name}_bare.svg > /dev/null
  png work/${name}_bare.svg figures/$name.png 700
}

# Étapes : Voronoï brut, Lloyd + fusion seule, réseau hexagonal + fusion/étirement
demo stage_raw      lattice=random lloydIterations=0  minEdgeRatio=0   minSides=3
demo stage_lloyd    lattice=random lloydIterations=10 minEdgeRatio=0.3 minSides=3
demo stage_hex      lattice=hex disorder=0.3 lloydIterations=5 minEdgeRatio=0.4 minSides=5

# Désordre du réseau hexagonal
for d in 0 0.3 0.5; do demo disorder_$d lattice=hex disorder=$d lloydIterations=5 minEdgeRatio=0.4 minSides=5; done

# Longueur de l'interface en pointe
for ti in 0 1 2; do demo tip_$ti lattice=hex disorder=0.3 lloydIterations=5 minEdgeRatio=0.4 minSides=5 tipInterface=$ti; done

# Exemple complet (géométrie de Frank-ouverture-fissure, paramètres fixes du tutoriel)
$NG tutorial_params.txt > work/tutorial.log
cp work/nodegen_input.txt figures/nodegen_input.txt
cp work/tutorial.log figures/console.txt
diff modele_input.txt work/input.txt > figures/deck_diff.txt || true
$NV work/nodefile.txt -bare -px 1000 -o work/example_bare.svg > /dev/null
png work/example_bare.svg figures/example.png 1000
$NV work/nodefile.txt -bare -nodes -box 0.0040 0.0048 0.0096 0.0104 -px 800 -o work/example_tip.svg > /dev/null
png work/example_tip.svg figures/example_tip.png 800
$NV work/nodefile.txt -bare -nodes -box 0.00319 0.00324 0.009975 0.010025 -px 800 -o work/example_mouth.svg > /dev/null
png work/example_mouth.svg figures/example_mouth.png 800
$NV work/nodefile.txt -px 1100 -o work/example_full.svg > /dev/null
png work/example_full.svg figures/example_legend.png 1100

# nodeFile de Frank (cellPrepro) pour comparaison
if [ -f ../../Frank-ouverture-fissure/nodefile.txt ]; then
  $NV ../../Frank-ouverture-fissure/nodefile.txt -bare -px 1000 -o work/frank_bare.svg > /dev/null
  png work/frank_bare.svg figures/frank.png 1000
fi
echo "figures régénérées dans $(pwd)/figures"
