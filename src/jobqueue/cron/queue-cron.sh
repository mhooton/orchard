#!/bin/bash
# Cron entry point for the orchard job queue on appct (user speculoos). See docs/job-queue.md.
#
#   queue-cron.sh [--shadow] watchdog | watch | lookback | report
#
# Runs one step inside orchard-server with docker exec, at most one copy at a time (flock), and
# appends its output to <queue_root>/logs/<step>.log. --shadow uses the shadow snapshot and the
# shadow database under mh_scratch, and never touches the container or the production tree.
#
#   watchdog  every 5 min: start a dispatcher if none holds the lock; SIGTERM one that stopped polling
#   watch     every 30 min through the afternoon: one ESO watcher poll
#   lookback  daily: the 30-night look-back
#   report    daily: write reports/report_<date>.md (live) or reports/shadow_<date>.md (shadow)

CONTAINER=orchard-server
COMPOSE=/appct/data/speculoos/orchard/docker/docker-compose.server.yml
HOST_DATA=/export/data/SPECULOOSPipeline      # is /data/SPECULOOSPipeline inside the container
SHADOW_HOME=/data/SPECULOOSPipeline/mh_scratch/appct-job-queue

if [ "$1" = "--shadow" ]; then
    shift
    MODE=shadow
    ROOT=$SHADOW_HOME/shadow
    EXEC_ENV=(-e "ORCHARD_QUEUE_CONFIG=$ROOT/config.json" -e "PYTHONPATH=$SHADOW_HOME/src:/opt/orchard/src")
else
    MODE=live
    ROOT=/data/SPECULOOSPipeline/queue
    EXEC_ENV=(-e "ORCHARD_QUEUE_ROOT=$ROOT")
fi
STEP=$1
case "$STEP" in
    watchdog|watch|lookback|report) ;;
    *) echo "usage: $0 [--shadow] watchdog|watch|lookback|report" >&2; exit 2 ;;
esac

HOST_ROOT=$HOST_DATA${ROOT#/data/SPECULOOSPipeline}
mkdir -p "$HOST_ROOT/logs" "$HOST_ROOT/reports"
exec 9>"/tmp/orchard-jobqueue-$MODE-$STEP.lock"
flock -n 9 || exit 0
exec >> "$HOST_ROOT/logs/$STEP.log" 2>&1

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
q() { docker exec "${EXEC_ENV[@]}" "$CONTAINER" python -m jobqueue "$@"; }

case "$STEP" in
    watchdog)
        if [ "$(docker inspect -f '{{.State.Running}}' "$CONTAINER" 2>/dev/null)" != "true" ]; then
            if [ "$MODE" = "shadow" ]; then
                echo "$(ts) $CONTAINER is not running; shadow mode leaves it alone"
                exit 0
            fi
            echo "$(ts) $CONTAINER is not running; docker-compose up -d"
            docker-compose -f "$COMPOSE" up -d || exit 1
        fi
        q watchdog
        rc=$?
        if [ $rc -eq 3 ]; then
            docker exec -d "${EXEC_ENV[@]}" "$CONTAINER" python -m jobqueue dispatcher --log
            echo "$(ts) started a $MODE dispatcher"
        elif [ $rc -ne 0 ] && [ $rc -ne 4 ]; then
            echo "$(ts) watchdog check failed (exit $rc)"
        fi
        ;;
    watch)
        q watch
        ;;
    lookback)
        q lookback
        ;;
    report)
        NAME=report
        [ "$MODE" = "shadow" ] && NAME=shadow
        q report --days 3 --out "$ROOT/reports/${NAME}_$(date -u +%Y%m%d).md" > /dev/null
        ;;
esac
