#!/usr/bin/env bash
# wrapper around mrbayes for bayesian inference tree building (Markov chain Monte Carlo analysis)
# consumes nexus produced by previous step, creates tree.nwk and stats.json

set -euo pipefail

# parse arguments
NEXUS=""
BURNINFRAC=""
OUT_TREE=""
OUT_STATS=""
SEED=""

while [[ $# -gt 0 ]]; do # loop while there are positional args left
  case "$1" in # look at first "--flag" in list
    --nexus)       NEXUS="$2"; shift 2;; # shift 2 = discard current "--flag value" and shifts $1 to next flag
    --burninfrac)  BURNINFRAC="$2"; shift 2;;
    --out-tree)    OUT_TREE="$2"; shift 2;;
    --out-stats)   OUT_STATS="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
  esac
done
# missing required args will cause an immediate code 2 exit here

# various specs
MB_NGEN=${MB_NGEN:-200000}
MB_SAMPLEFREQ=${MB_SAMPLEFREQ:-200}
MB_PRINTFREQ=${MB_PRINTFREQ:-1000}
MB_DIAGNFREQ=${MB_DIAGNFREQ:-5000}
MB_NCHAINS=${MB_NCHAINS:-4}
MB_TEMP=${MB_TEMP:-0.1}
MB_STOPRULE=${MB_STOPRULE:-yes}
MB_STOPVAL=${MB_STOPVAL:-0.01}

# small checks
[[ -n "$NEXUS" && -n "$BURNINFRAC" && -n "$OUT_TREE" && -n "$OUT_STATS" ]] || {
  echo "[run_mrbayes.sh] missing required arguments" >&2
  exit 2
}
[[ -f "$NEXUS" ]] || { echo "[run_mrbayes.sh] nexus not found: $NEXUS" >&2; exit 2; }

# pick MrBayes executable, prefer mb
if command -v mb >/dev/null 2>&1; then
  MBEXE=mb
elif command -v mb-mpi >/dev/null 2>&1; then
  MBEXE=mb-mpi
else
  echo "[run_mrbayes.sh] MrBayes executable not found (looked for 'mb' or 'mb-mpi')" >&2
  exit 2
fi

# output prefix for the files mrbayes will generate
PREFIX="${OUT_TREE%.*}"
OUT_DIR="$(dirname "$OUT_TREE")"
mkdir -p "$OUT_DIR" "$(dirname "$OUT_STATS")"

# remove previous outputs with same prefix to avoid interactive prompts
rm -f "${PREFIX}".run*.p "${PREFIX}".run*.t \
      "${PREFIX}".ckp "${PREFIX}".ckp~ "${PREFIX}".mcmc \
      "${PREFIX}".trprobs "${PREFIX}".tstat "${PREFIX}".vstat \
      "${PREFIX}".parts "${PREFIX}".con.tre 2>/dev/null || true

# run MrBayes in batch
start_ts=$(date +%s)

set +e # capture mrbayes exit code for troubleshooting

# batch commands sent to MrBayes executable:
# - execute nexus
# - set filenames (prefix)
# - optionally sets seed
# - applies MCMC params and burninfrac
# - summarize into consensus tree
# - quit 

"$MBEXE" <<MB
execute $NEXUS;
mcmcp filename=$PREFIX;
$( [[ -n "$SEED" ]] && echo "set seed=$SEED swapseed=$SEED;" )
mcmcp ngen=$MB_NGEN samplefreq=$MB_SAMPLEFREQ printfreq=$MB_PRINTFREQ diagnfreq=$MB_DIAGNFREQ nchains=$MB_NCHAINS temp=$MB_TEMP stoprule=$MB_STOPRULE stopval=$MB_STOPVAL;
mcmcp burninfrac=$BURNINFRAC;
mcmc;
sumt burninfrac=$BURNINFRAC contype=allcompat;
quit;
MB


rc=$? # exit code of the most recent command

set -e # re-enable fast failing

#fail explicitly if mrbayes failed

if [[ $rc -ne 0 ]]; then
  echo "[run_mrbayes.sh] MrBayes failed with exit code $rc" >&2
  exit $rc
fi

end_ts=$(date +%s)

elapsed=$(( end_ts - start_ts ))

################### collect outputs ###############################

# MrBayes writes a consensus tree to <prefix>.con.tre by default
CON_TRE="${PREFIX}.con.tre"
if [[ ! -s "$CON_TRE" ]]; then
  echo "[run_mrbayes.sh] expected consensus tree not found: $CON_TRE" >&2
  exit 2
fi

# copy final Newick to requested location
cp -f "$CON_TRE" "$OUT_TREE"

# write stats json 
cat > "$OUT_STATS" <<JSON
{
  "method": "BI",
  "burninfrac": ${BURNINFRAC},
  "ngen": ${MB_NGEN},
  "samplefreq": ${MB_SAMPLEFREQ},
  "printfreq": ${MB_PRINTFREQ},
  "diagnfreq": ${MB_DIAGNFREQ},
  "nchains": ${MB_NCHAINS},
  "temp": ${MB_TEMP},
  "stoprule": "$(printf '%s' "$MB_STOPRULE")",
  "stopval": ${MB_STOPVAL},
  "elapsed_sec": ${elapsed}
}
JSON
