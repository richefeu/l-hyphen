#!/usr/bin/env bash
# =====================================================================================================
# run_cases.sh : lance les cas d'une étude paramétrique l-hyphen en parallèle, avec reprise.
#
# Les cas sont les répertoires listés dans <cases>/cases.txt (produit par generate_cases), chacun contenant
# un input.txt. Chaque cas est lancé dans son répertoire (« run input.txt > log.txt ») ; son état est noté
# dans <cas>/status.txt :
#
#   done     <date> <durée> s  <fin>     calcul terminé (fin = « événement » ou « TFINAL »)
#   failed   <date> code <rc>  <raison>  échec : code de retour non nul, durée dépassée, aucune cellule (nodeFile
#                                        introuvable : run rend alors 0), nan dans les résultats
#   running  <hôte> <pid> <date>         calcul en cours (ou interrompu, si le processus n'existe plus)
#
# Relancer le script reprend l'étude : les cas « done » sont sautés, les cas en échec ou interrompus sont
# relancés depuis le début (sorties précédentes effacées), dans la limite de -m tentatives. Un cas en cours
# sur une autre machine (répertoire partagé) n'est pas touché. Ctrl-C arrête proprement tous les calculs ;
# ils seront relancés à la reprise.
#
# Compatible bash 3.2 (macOS) et Linux.
# =====================================================================================================

usage() {
  cat <<EOF
Usage : $0 [options]
  -d DIR   répertoire des cas (défaut : cases)
  -j N     nombre de calculs simultanés (défaut : nombre de cœurs)
  -r RUN   exécutable l-hyphen (défaut : \$LHYPHEN_RUN, sinon ../run, sinon « run » du PATH)
  -t SEC   durée maximale d'un calcul en secondes (défaut : aucune)
  -m N     nombre maximal de tentatives par cas (défaut : 3)
  -s       affiche seulement l'état de l'étude
  -h       aide
EOF
}

script_dir=$(cd "$(dirname "$0")" && pwd)
cases_dir="cases"
jobs_max=""
run_bin="${LHYPHEN_RUN:-}"
time_max=0
attempts_max=3
status_only=0

while getopts "d:j:r:t:m:sh" opt; do
  case $opt in
    d) cases_dir=$OPTARG ;;
    j) jobs_max=$OPTARG ;;
    r) run_bin=$OPTARG ;;
    t) time_max=$OPTARG ;;
    m) attempts_max=$OPTARG ;;
    s) status_only=1 ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

if [ ! -d "$cases_dir" ]; then
  echo "run_cases : répertoire des cas introuvable : $cases_dir (lancer d'abord ./generate_cases)" >&2
  exit 1
fi
cases_dir=$(cd "$cases_dir" && pwd)
host=$(hostname)

# --- liste des cas : cases.txt (2e colonne), sinon tous les sous-répertoires contenant un input.txt
case_list() {
  if [ -f "$cases_dir/cases.txt" ]; then
    grep -v '^[[:space:]]*#' "$cases_dir/cases.txt" | awk 'NF >= 2 { print $2 }'
  else
    for d in "$cases_dir"/*/; do [ -f "$d/input.txt" ] && basename "$d"; done
  fi
}

# --- état d'un cas : done, failed, running (processus vivant ici), remote (en cours sur une autre machine),
#     stale (en cours mais processus disparu), todo
case_state() {
  local st="$cases_dir/$1/status.txt" w h p
  [ -f "$st" ] || { echo todo; return; }
  read -r w h p _ < "$st"
  case $w in
    done) echo done ;;
    failed) echo failed ;;
    running)
      if [ "$h" != "$host" ]; then echo remote
      elif kill -0 "$p" 2>/dev/null; then echo running
      else echo stale
      fi ;;
    *) echo todo ;;
  esac
}

attempts_of() {
  local f="$cases_dir/$1/attempts.txt"
  if [ -f "$f" ]; then cat "$f"; else echo 0; fi
}

print_status() {
  local n_done=0 n_failed=0 n_running=0 n_remote=0 n_stale=0 n_todo=0 c s
  for c in $(case_list); do
    s=$(case_state "$c")
    case $s in
      done) n_done=$((n_done + 1)) ;;
      failed) n_failed=$((n_failed + 1)); echo "  échec      : $c ($(cut -d' ' -f2- "$cases_dir/$c/status.txt"), tentatives $(attempts_of "$c"))" ;;
      running) n_running=$((n_running + 1)); echo "  en cours   : $c" ;;
      remote) n_remote=$((n_remote + 1)); echo "  en cours   : $c (sur $(cut -d' ' -f2 "$cases_dir/$c/status.txt"))" ;;
      stale) n_stale=$((n_stale + 1)); echo "  interrompu : $c" ;;
      *) n_todo=$((n_todo + 1)) ;;
    esac
  done
  echo "État : $n_done terminés, $n_failed en échec, $((n_running + n_remote)) en cours, $n_stale interrompus, $n_todo à faire"
}

if [ "$status_only" = 1 ]; then
  print_status
  exit 0
fi

# --- exécutable l-hyphen
if [ -z "$run_bin" ]; then
  if [ -x "$script_dir/../run" ]; then run_bin="$script_dir/../run"
  else run_bin=$(command -v run || true)
  fi
fi
if [ -z "$run_bin" ] || [ ! -x "$run_bin" ]; then
  echo "run_cases : exécutable l-hyphen introuvable (option -r ou variable LHYPHEN_RUN)" >&2
  exit 1
fi
run_bin=$(cd "$(dirname "$run_bin")" && pwd)/$(basename "$run_bin")

# --- nombre de calculs simultanés
if [ -z "$jobs_max" ]; then
  jobs_max=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
fi

# --- un cas : nettoyage, calcul, contrôle du résultat
run_case() {
  local c=$1 dir="$cases_dir/$1" n start rc rpid wd reason end dur
  n=$(( $(attempts_of "$c") + 1 ))
  echo "$n" > "$dir/attempts.txt"
  ( cd "$dir" && "$run_bin" -c > /dev/null 2>&1 )        # efface conf*, sample*.svg, breakEvol.txt, ...
  rm -f "$dir/top.txt" "$dir/bottom.txt" "$dir/log.txt"
  start=$(date +%s)

  cd "$dir" || return
  "$run_bin" input.txt > log.txt 2>&1 &
  rpid=$!
  echo "running $host $rpid $(date '+%Y-%m-%dT%H:%M:%S')" > status.txt
  wd=""
  if [ "$time_max" -gt 0 ] 2>/dev/null; then
    ( sleep "$time_max"; kill "$rpid" 2>/dev/null ) &
    wd=$!
  fi
  wait "$rpid"
  rc=$?
  if [ -n "$wd" ]; then
    pkill -P "$wd" 2>/dev/null   # le sleep du chien de garde
    kill "$wd" 2>/dev/null
  fi
  dur=$(( $(date +%s) - start ))

  reason=""
  if [ "$rc" -ne 0 ]; then
    if [ "$time_max" -gt 0 ] 2>/dev/null && [ "$dur" -ge "$time_max" ]; then reason="durée maximale dépassée"
    else reason="code de retour non nul"
    fi
  elif [ ! -f diagnostic.txt ]; then
    reason="pas de diagnostic.txt (erreur de lecture de l'input ?)"
  elif ! grep -Eq 'Cellules +: +[1-9]' diagnostic.txt; then
    reason="aucune cellule (nodeFile introuvable ?)"
  elif cat top.txt bottom.txt breakEvol.txt 2>/dev/null | grep -qi 'nan'; then
    reason="nan dans les résultats (divergence)"
  fi

  if [ -z "$reason" ]; then
    if grep -q "end of simulation" log.txt; then end="événement"; else end="TFINAL (sans événement d'arrêt)"; fi
    echo "done $(date '+%Y-%m-%dT%H:%M:%S') $dur s $end" > status.txt
    echo "[$(date '+%H:%M:%S')] terminé  : $c ($dur s, fin : $end)"
  else
    echo "failed $(date '+%Y-%m-%dT%H:%M:%S') code $rc $reason" > status.txt
    echo "[$(date '+%H:%M:%S')] ÉCHEC    : $c (tentative $n, $reason ; voir $c/log.txt)"
  fi
}

# --- Ctrl-C / arrêt : on arrête les calculs lancés par ce script (et seulement eux : pas de « kill 0 », qui
#     toucherait aussi le shell appelant) ; les cas restent « running » avec un processus disparu, ils seront
#     donc relancés à la reprise
stop_all() {
  trap - INT TERM
  echo
  echo "Interruption : arrêt des calculs en cours (ils seront relancés à la reprise)."
  local j
  for j in $(jobs -p); do
    pkill -TERM -P "$j" 2>/dev/null   # run et chien de garde lancés par le cas
    kill -TERM "$j" 2>/dev/null
  done
  wait
  exit 130
}
trap stop_all INT TERM

echo "run_cases : $(case_list | wc -l | tr -d ' ') cas dans $cases_dir, $jobs_max calculs simultanés, exécutable $run_bin"
t0=$(date +%s)
n_launched=0
for c in $(case_list); do
  if [ ! -f "$cases_dir/$c/input.txt" ]; then
    echo "  ignoré     : $c (pas d'input.txt)"
    continue
  fi
  s=$(case_state "$c")
  case $s in
    done|running|remote) continue ;;
    failed|stale)
      if [ "$(attempts_of "$c")" -ge "$attempts_max" ]; then
        echo "  abandonné  : $c ($(attempts_of "$c") tentatives, voir $c/log.txt)"
        continue
      fi ;;
  esac
  while [ "$(jobs -rp | wc -l)" -ge "$jobs_max" ]; do sleep 1; done
  echo "[$(date '+%H:%M:%S')] lancement : $c$( [ "$s" = todo ] || echo " (reprise : $s)")"
  run_case "$c" &
  n_launched=$((n_launched + 1))
done
wait
trap - INT TERM

echo "run_cases : $n_launched cas lancés en $(( $(date +%s) - t0 )) s"
print_status
