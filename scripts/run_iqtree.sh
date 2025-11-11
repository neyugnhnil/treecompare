#!/usr/bin/env bash
# wrapper around IQ-TREE for Maximum Likelihood tree building

set -euo pipefail

################ parse arguments from dispatcher.py ######################
FASTA=""
DATATYPE=""
MODEL=""
STARTER=""
OUT_TREE=""
OUT_STATS=""
SEED=""

while [[ $# -gt 0 ]]; do # loop while there are positional args left
  case "$1" in # look at first "--flag" in list
    --fasta) FASTA="$2"; shift 2;; # shift 2 = discard current "--flag value" and shifts $1 to next flag
    --datatype) DATATYPE="$2"; shift 2;;
    --model) MODEL="$2"; shift 2;;
    --starter) STARTER="$2"; shift 2;;
    --out-tree) OUT_TREE="$2"; shift 2;;
    --out-stats) OUT_STATS="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
  esac
done
# missing required args will cause an immediate code 2 exit here

########## quick recheck of validation done upstream #####################
[[ -n "$FASTA" && -n "$DATATYPE" && -n "$MODEL" && -n "$STARTER" && -n "$OUT_TREE" && -n "$OUT_STATS" ]] || {
  echo "[run_iqtree.sh] missing required arguments" >&2
  exit 2
}

############# prefers iqtree2 ############################################
if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE=iqtree2
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE=iqtree
else
  echo "[run_iqtree.sh] iqtree executable not found (looked for iqtree2 or iqtree)" >&2
  exit 2
fi

############ misc. IQ-TREE syntax #########################################

# map datatype to IQ-TREE -st values 
case "$DATATYPE" in
  dna|DNA) STYPE="DNA" ;;
  aa|AA) STYPE="AA" ;;
  *) echo "[run_iqtree.sh] datatype must be dna or aa" >&2; exit 2 ;; # just in case
esac

# IQ-TREE needs a prefix for outputs not a path
PREFIX="${OUT_TREE%.*}"

# ensure output dir exists
mkdir -p "$(dirname "$OUT_TREE")" "$(dirname "$OUT_STATS")"

# starter: `-t <file>` for a starting tree; `-t RANDOM` for random start.
# as array to be expanded
if [[ "${STARTER,,}" == "random" ]]; then
  START_OPT=(-t RANDOM)
else
  START_OPT=(-t "$STARTER")
fi

# optional seed if provided
# as array to be expanded (expands to nothing if empty)
SEED_OPT=() 
if [[ -n "$SEED" ]]; then
  SEED_OPT=(-seed "$SEED")
fi

########## assemble command that runs IQ-TREE ############################
start_ts=$(date +%s)
set +e
"$IQTREE" \
  -s "$FASTA" -m "$MODEL" -st "$STYPE" "${START_OPT[@]}" \
  -nt AUTO -pre "$PREFIX" -quiet \
  "${SEED_OPT[@]}"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "[run_iqtree.sh] IQ-TREE failed with exit code $rc" >&2
  exit $rc
fi
end_ts=$(date +%s)
elapsed=$(( end_ts - start_ts ))

################### collect outputs ###############################
# write ${PREFIX}.treefile
TREEFILE="${PREFIX}.treefile"
if [[ ! -s "$TREEFILE" ]]; then
  echo "[run_iqtree.sh] expected treefile not found: $TREEFILE" >&2
  exit 2
fi

# copy tree to requested path
cp -f "$TREEFILE" "$OUT_TREE"

# collect stats from .iqtree text report if available
IQREP="${PREFIX}.iqtree"
best_ll=""
model_used=""
if [[ -f "$IQREP" ]]; then
  # best score found
  best_ll=$(grep -E "BEST SCORE FOUND" "$IQREP" | awk -F':' '{print $2}' | xargs || true)
  # model of substitution
  model_used=$(grep -E "Model of substitution" "$IQREP" | awk -F':' '{print $2}' | xargs || true)
fi

# stats JSON 
cat > "$OUT_STATS" <<JSON
{
  "method": "ML",
  "datatype": "${STYPE}",
  "model": "$(printf '%s' "${model_used:-$MODEL}")",
  "starter": "$( [[ "${STARTER,,}" == "random" ]] && echo random || echo provided )",
  "elapsed_sec": ${elapsed},
  "iqtree_exe": "${IQTREE}",
  "prefix": "$(printf '%s' "$PREFIX")",
  "best_loglik": "$(printf '%s' "${best_ll:-}")",
  "outputs": {
    "treefile": "$(printf '%s' "$TREEFILE")",
    "iqtree_report": "$( [[ -f "$IQREP" ]] && printf '%s' "$IQREP" || echo "" )"
  }
}
JSON
