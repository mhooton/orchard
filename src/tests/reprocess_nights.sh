#!/bin/bash
# reprocess_nights.sh — run a list of nights of one target through the
# pipeline, sequentially, inside the pipeline container.
#
#   reprocess_nights.sh NIGHTS_FILE TARGET [ZLP flags...]
#
# NIGHTS_FILE has one "TELESCOPE DATE" per line (a leading '#' comments a
# line out).  Every night is run with the same flags; a failing night is
# logged and the loop continues.  Per-night stdout goes to
# $LOGDIR/<TEL>_<DATE>.log, a one-line summary per night to $LOGDIR/summary.txt.
#
# Environment (all optional):
#   ORCHARD_PATH       pipeline source dir (default: this script's ../ )
#   BASEDIR            data root inside the container (default /data/SPECULOOSPipeline)
#   GAIADATABASEPATH   Gaia database
#   TARGET_LIST        40 pc target list override
#   N_CORES            cores per night (default 20)
#   LOGDIR             where to write logs (default $ORCHARD_PATH/../reprocess_logs)
#
# Example, launched detached from the container host:
#   docker exec -d orchard-server bash -c "cd /gaia_database/orchard-missing-gaia/src && \
#     GAIADATABASEPATH=/gaia_database/gaia_dr3_unified_16jcut_dr3only.db \
#     TARGET_LIST=/gaia_database/orchard-missing-gaia/ml_40pc.txt \
#     bash tests/reprocess_nights.sh ../nights_Io.txt Sp1056+0700 --no_T12 --force-platesolve"
set -u
NIGHTS=$1; TARGET=$2; shift 2; FLAGS="$*"
HERE=$(cd "$(dirname "$0")" && pwd)
export ORCHARD_PATH=${ORCHARD_PATH:-$(dirname "$HERE")}
export PYTHONPATH=$ORCHARD_PATH
export N_CORES=${N_CORES:-20}
BASEDIR=${BASEDIR:-/data/SPECULOOSPipeline}
LOGDIR=${LOGDIR:-$(dirname "$ORCHARD_PATH")/reprocess_logs}
mkdir -p "$LOGDIR"
cd "$ORCHARD_PATH"
n_ok=0; n_fail=0
while read -r TEL DATE _; do
    [ -z "${TEL:-}" ] && continue
    case "$TEL" in \#*) continue ;; esac
    LOG=$LOGDIR/${TEL}_${DATE}.log
    t0=$(date +%s)
    echo "[$(date -u +%FT%TZ)] START $TEL $DATE -> $LOG" | tee -a "$LOGDIR/summary.txt"
    if ./main/ZLP_pipeline.sh $FLAGS 1 "$BASEDIR" "$DATE" 8 2 "$TEL" "$TARGET" > "$LOG" 2>&1; then
        status=OK; n_ok=$((n_ok+1))
    else
        status="FAILED(exit $?)"; n_fail=$((n_fail+1))
    fi
    dt=$(( $(date +%s) - t0 ))
    nmcmc=$(ls "$BASEDIR/PipelineOutput/v3/$TEL/output/$DATE/$TARGET/lightcurves/"*_MCMC 2>/dev/null | wc -l)
    echo "[$(date -u +%FT%TZ)] END   $TEL $DATE $status ${dt}s mcmc_files=$nmcmc" | tee -a "$LOGDIR/summary.txt"
done < "$NIGHTS"
echo "[$(date -u +%FT%TZ)] DONE ok=$n_ok failed=$n_fail" | tee -a "$LOGDIR/summary.txt"
