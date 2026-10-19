#!/bin/bash
# Régénère les figures du tutoriel nodegen (nécessite inkscape).
# Usage : cd doc && ./make_figures.sh
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
demo stage_raw      lattice=random lloydIterations=0  minEdgeFactor=0   minSides=3
demo stage_lloyd    lattice=random lloydIterations=10 minEdgeFactor=0.19 minSides=3
demo stage_hex      lattice=hex disorder=0.3 lloydIterations=5 minEdgeFactor=0.25 minSides=5

# Désordre du réseau hexagonal
for d in 0 0.3 0.5; do demo disorder_$d lattice=hex disorder=$d lloydIterations=5 minEdgeFactor=0.25 minSides=5; done

# Longueur de l'interface en pointe
for ti in 0 1 2; do demo tip_$ti lattice=hex disorder=0.3 lloydIterations=5 minEdgeFactor=0.25 minSides=5 tipInterface=$ti; done

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
# ---------------------------------------------------------------------------------------------------------
# Tissus végétaux : bundles et attraction de Mahalanobis
# ---------------------------------------------------------------------------------------------------------
tissu() { # nom puis paires clé=valeur ; domaine de 8 mm, ~1300 cellules, 5 % de bundles
  local name=$1; shift
  cat > work/$name.params <<EOF
Lx 0.008
Ly 0.008
cellSize 2.5e-4
lattice random
lloydIterations 30
bundleFraction 0.05
seed 2
output work/$name.txt
svg work/$name.svg
input work/${name}_input.txt
bundleOutput work/${name}_bundles.txt
EOF
  for kv in "$@"; do echo "${kv%%=*} ${kv#*=}" >> work/$name.params; done
  $NG work/$name.params > work/$name.log
  $NV work/$name.txt -bare -w 2e-6 -px 700 -fill work/${name}_bundles.txt -o work/${name}_bare.svg > /dev/null
  png work/${name}_bare.svg figures/$name.png 700
}
for a in 0 0.05 0.2; do tissu tissu_alpha_$a alpha=$a; done
for K in 1 100 1000; do tissu tissu_K_$K alpha=0.05 anisotropy=$K; done
tissu tissu_tomo lattice=square disorder=0 lloydIterations=5 alpha=0.02 anisotropy=10000
tissu tissu_fixed alpha=0.05 bundleMotion=fixed
tissu tissu_attract alpha=0.05 bundleMotion=attract

# Balayages (domaine de 1 cm, ~2000 cellules) : une ligne de métriques par calcul (mot-clé morphology)
sweep() { # fichier de métriques puis paires clé=valeur
  local out=$1; shift
  cat > work/sweep.params <<EOF
Lx 0.01
Ly 0.01
cellSize 2.5e-4
lattice random
lloydIterations 30
bundleFraction 0.05
seed 1
output work/sweep.txt
svg work/sweep.svg
input work/sweep_input.txt
bundleOutput work/sweep_bundles.txt
morphology $out
EOF
  for kv in "$@"; do echo "${kv%%=*} ${kv#*=}" >> work/sweep.params; done
  $NG work/sweep.params > /dev/null
}
rm -f work/sweep_*.dat
for a in 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2; do
  sweep work/sweep_alpha.dat alpha=$a
done
for K in 1 2 5 10 20 50 100 200 500 1000; do sweep work/sweep_K.dat alpha=0.05 anisotropy=$K; done
for n in 1 2 5 10 20 30 50 100 200; do
  sweep work/sweep_it_K1.dat alpha=0.05 lloydIterations=$n
  sweep work/sweep_it_K1000.dat alpha=0.05 anisotropy=1000 lloydIterations=$n
done
cp work/sweep_*.dat figures/
cp work/tissu_alpha_0.05.log figures/console_tissu.txt

echo "figures régénérées dans $(pwd)/figures"
