# The appct job queue

A priority queue for the orchard pipeline on appct (`src/jobqueue/`), an ESO watcher that queues each night as
soon as its frames are at ESO, and a daily look-back over the last 30 nights. It is decision Q1 of the portal
completeness plan, with Q3 and Q4 (nightly processing through the queue, download separate from processing) as
the first use. Backlog reruns (class P2) come later, once the 16-bit frames, a per-run catcache and the K8 fixes
have landed.

Everything is standard-library Python 3.9 that runs inside `orchard-server` on appct. The ESO side reuses
`download.request_eso.ESODownloader` and, for the nightly download, `download/SSO_download.py` itself.

## Why

The completeness audit of 8–9 October 2026 found three things this addresses:

- **Latency.** A night's frames reach ESO 3–4 h after it ends (13:30–15:00 UTC), but the download only starts at
  a fixed cron slot (19:00–22:00 BST). The median from last frame to portal was 13.7 h, of which about 4–6 h is
  waiting for the slot.
- **One attempt, never revisited.** The cron fetches yesterday only. Frames that reach ESO late, or nights the
  cron missed, stay missing; a download error also skips that night's processing (`&&`). Io has been down since
  29 August and still retries 24 times a night.
- **An idle machine with a shared bottleneck.** appct averaged 0.8 % CPU, but its two spinning arrays saturate
  with five bulk jobs, and nothing coordinates who uses them.

## How it works

### Classes

| Class | What | Source |
|---|---|---|
| P0 | last night: its download, then its pipeline run | the ESO watcher |
| P1 | a recent night that gained frames, or has no light curves | the look-back, and fetch jobs |
| P2 | backlog work: old nights to download and/or reprocess | `q add` (`--download`, `--from-file`) |
| P3 | background: 16-bit sweep, PWV-only replays, the nightly ledger (later) | `q add` |

### Admission

Each job declares its cores and whether it is disk-heavy. Every poll (45 s) the dispatcher starts a queued job
only if all of these hold:

- declared cores of running jobs plus its own fit in **112** of the 128 threads; from 11:00 to 01:00 UTC,
  **40** of those are kept for P0, so P1–P3 get at most 72;
- the 5-minute load average, plus the cores started earlier in the same poll, is under **100**;
- `MemAvailable` is over **64 GB**;
- both arrays (`sda` = PipelineOutput, `sdb` = data/raw) are under **70 %** utilisation, measured from
  `/proc/diskstats` over the last poll interval. P0 skips this check. At most one disk-heavy job starts per
  poll, because `%util` only shows a new job's I/O at the next poll;
- its locks are free (below);
- a P1–P3 job on a telescope stays off it if that telescope's nightly job is due (13:00–17:00 UTC for SSO,
  10:00–13:00 UTC for Artemis) before the job would finish (`p0_guard`).

### Strict priority, never preemption

Jobs are considered by class, then priority within the class, then age. When a job cannot start, the resources
it is waiting for are reserved for its class: no lower-class job may take them. A lower-class job may still use
something else (a P1 job on Europa runs while a P0 job waits for the Callisto lock), and a same-class job may use
what is left (a P0 download is not held up by a P0 pipeline job waiting for cores).

A running job is never killed or paused to make room: the plate solver relies on a 60 s `SIGALRM`, and the
pipeline holds files open. The P0 reserve and the P0 guard keep the wait short instead. The only kills are a job
past three times its estimate (floor 120 min) and an explicit `q cancel --kill`.

### Locks

| Lock | Held by | Why |
|---|---|---|
| `tel:<TEL>` | pipeline jobs | T8 and T9 run `rm -rf ${DATDIR}/catcache`, one directory per telescope. Off (`telescope_lock: false`) once a per-run catcache lands (K7). |
| `night:<TEL>:<NIGHT>` | every job on that night | raw and output directories of one night |
| `target:<NAME>` | pipeline jobs, per target | `PipelineOutput/v2/StackImages` is shared across telescopes and nights |
| `sem:eso` | download and fetch jobs | at most `semaphores.eso` (2) at once |
| `sem:eso_bulk` | P2 and P3 fetch jobs as well | at most `semaphores.eso_bulk` (1), so bulk downloads never hold both ESO slots and the nightly download always finds one free |

Pipeline and download processes the queue did not start (the cron, a person, a replay) hold the same locks: the
dispatcher reads them from `/proc/*/cmdline` inside the container. Processes in other containers
(`orchard-server-spirit`, `orchard-server-autosave`) are invisible to it.

Targets for locking come from ESO's `object` column (SSO) or the raw file names (Artemis), upper-cased with
spaces as `--`, as `createlists.py` names them.

### Jobs, logs, retries

Each job runs through `jobqueue.runner` in its own session, so the job and everything it starts share one
process group. Output goes to `<queue_root>/logs/jobs/<YYYY-MM>/<id>-<kind>-<tel>-<night>-a<attempt>.log`, which
unlike `PipelineOutput/v3/<TEL>/logs/<NIGHT>_1_v3.log` is never overwritten by a later run. The runner writes
`<log>.status.json` with the exit code when the command ends.

Every job records its start, end, exit code and the deployed code version: the first word of `src/VERSION` if
it exists, otherwise `tree-<hash>` of the deployed `.py` and `.sh` files (the production tree is a working copy,
not a git checkout).

When a job fails, its **failure signature** is the last exception line of its log (else the last error line),
with paths, numbers and target names normalised. It is **transient** if the signature or the last 30 lines
match a network, full-disk, I/O or memory pattern, or the exit code is 75, 137 or 143; transient failures are
retried after an hour, at most twice. Anything else is **deterministic**: the job stops as failed, and after a
fix `q requeue --signature "<text>"` puts every job with that signature back.

### Restarts and the watchdog

The dispatcher holds an `flock` on `<queue_root>/dispatcher.lock` for its whole life, so only one runs. Every
5 minutes cron runs `queue-cron.sh watchdog`, which starts one (`docker exec -d`) if nobody holds the lock, and
sends `SIGTERM` to one that has not polled for 10 minutes.

Jobs survive a dispatcher restart: the new dispatcher finds them by PID (checked against the process start time,
so a reused PID is not mistaken for the job) and collects their exit status from the status file. A job whose
process vanished without a status file (container restart) is requeued once; the second time it fails as
`vanished`. If only the runner dies (the OOM killer, say) while its command carries on, the job stays running
until the whole process group has gone, so it is never started twice; it is no longer signalled, because its
group can no longer be verified.

### The database

SQLite in WAL mode at `<queue_root>/queue.sqlite`, on appct's local disk (`/export/data/SPECULOOSPipeline/queue`
= `/data/SPECULOOSPipeline/queue` in the container, XFS on `sdb1`, not the `sync`-mounted PipelineOutput).
appcs sees that disk over NFS, where SQLite locking is unreliable, so the store refuses to open on a network
filesystem: always use it through `docker exec orchard-server` on appct.

Tables: `jobs`, `events` (every transition, waiting reason, note and watcher decision), `watch`, `lookback`, `kv`
(pauses, drain, dispatcher heartbeat).

## The ESO watcher

`queue-cron.sh watch`, every 30 minutes from 11:00 to 23:30 local time. For last night and the night before, it
asks ESO for each SSO telescope's survey frames (`prog_id` 60.A-9009(A–D), `dp_id LIKE 'SPECU%'`, night D =
[D 15:00, D+1 15:00) UTC) through the pipeline's authenticated client. The credentials come from
`src/reporting/.env` and are never printed: ESODownloader's own messages can contain the token URL with the
password, so they are captured and scrubbed.

A night is **ready** when:

1. ESO's count reaches the count in `Observations/<TEL>/transfer_log.txt` (for Callisto SPIRIT both count
   datacubes); or
2. with no transfer-log entry, the count and ESO's latest `last_mod_date` are unchanged on two polls at least
   25 minutes apart, and nothing changed at ESO for 20 minutes; or
3. below the transfer-log count, nothing has changed for 3 hours (it queues what is there; the look-back fetches
   the rest).

Then it queues two P0 jobs: a **download** (`python -m jobqueue jobs download`) and the **pipeline** run after
it, exactly the cron's command (`./main/ZLP_pipeline.sh --force-platesolve 1 /data/SPECULOOSPipeline <NIGHT> 8 2
<TEL>`). The download job:

- runs `download/SSO_download.py` for that one night with `--max-retries 1`, logging to
  `ESO_logs/<NIGHT><TEL>.log` and sending its usual summary email, as the cron does;
- refuses any night more than 5 days old, so it can never be looped over old nights;
- first reads the headers of the frames already on disk and, if `transformation_check()` would delete any of
  them (Callisto all dates, Ganymede after 2025-03-03, frames without `ASTRAROT`/`ASTRAMIR`), does not run
  `SSO_download.py` at all but fetches the missing frames add-only;
- exits 0 whenever frames are on disk, so the pipeline runs even if a few are missing; exits 75 (retry in an hour)
  only if nothing landed.

Other cases:

- **No frames** at ESO and no transfer-log entry by 21:00 UTC the next day: one quiet note in the event log, no
  download, no email, no pipeline run (Io since 29 August). Frames that turn up later are still queued.
- **Artemis** does not go through ESO. Its night is ready when
  `Observations/Artemis/images/<NIGHT>/Data_Download.txt` (written at the
  end of each transfer, about 10:45–11:45 UTC) lists only files that are present, or when the directory has
  stopped changing for two polls; at 22:00 UTC it is queued with whatever is there. No download job.
- A night that already has a pipeline log (the cron or a person ran it) is marked `processed` and left to the
  look-back, so the switch-over never reprocesses the night before last.
- If ESO cannot be reached the watcher changes nothing; the next poll retries.

## The look-back

`queue-cron.sh lookback`, daily at 21:15 local. One ESO query per telescope lists every frame of the last
30 nights (last night excluded, it belongs to the watcher), and each night is compared with its directory:

- An ANDOR row is present when `<dp_id>.fits` (or `.fits.fz`, `.fts`, …) exists. A SPIRIT datacube is unpacked
  on download into frames named by their own timestamps, so a cube row (`origfile` `SPECU4.<YYYYMMDDTHHMMSS>_…`;
  `det_ndit` is empty) counts as present when any frame on disk falls between its start and the next cube's
  start. On Callisto 5–7 October every unpacked frame landed in a cube this way.
- **Frames missing** → a P1 **fetch** job (add-only, below). If it adds science frames it queues a P1 rerun of the
  night itself. Rows it downloaded but that added nothing are remembered and not asked for again.
- **Science frames on disk but no light curve in v2 or v3** → a P1 rerun of the night, at most twice ever.
- **Light curves only in v3** → reported (a T12 promotion question, never automatic).
- Nights with a queued or running job are left alone. At most 3 fetches per night, 3,000 rows per night and
  10,000 rows per run.

### The add-only fetch

`SSO_download.py` must never be looped over old nights: `transformation_check()` deletes every pre-Astra Callisto
frame before re-downloading it; it counts only `.fits`/`.fts` files as present, so it fetches the 691 `.fits.fz`
nights and every unpacked SPIRIT night all over again (and unpacking with `overwrite=True` wipes the WCS T7
wrote into those raw frames); it duplicates the 2018 ACP-named frames; and it creates empty night directories,
rewrites `download_log.csv` and emails every night. `jobqueue/fetch.py` is the one download path for anything
that is not last night, used by the look-back, by `q add --download` and by the backfill. It only adds:

- **what is missing** is decided by `nights.diff_eso_disk`: an ESO row is present if a file named after its
  dp_id exists with any FITS suffix (so `.fits.fz` counts), if it is a SPIRIT datacube and frames unpacked from
  it are on disk (see the look-back), or if a frame under another name, in the night directory or its
  `AutoFlat/` and `Calibration/` subdirectories, has the same `DATE-OBS`. ESO's dp_id is `SPECU<n>.<DATE-OBS>`:
  ACP's `Sp1609-3431-S001-R001-C001-I+z.fts` on Io 2018-09-05 has `DATE-OBS = '2018-09-05T23:11:08.710'` and is
  `SPECU1.2018-09-05T23:11:08.710` at ESO, with identical pixels (ESO stores the 2018 frames as float32, so
  their copies take twice the disk);
- each file is downloaded under a temporary name into `<queue_root>/staging/` (same filesystem as the archive),
  at most `fetch.workers` (3) at a time, each with `fetch.attempts` (4) tries and exponential backoff from 30 s;
  after 8 network failures in a row the job stops and exits 75, so the queue retries it an hour later. gzip
  and compress data are decompressed as `request_eso.py` does;
- each file is **verified** before anything else happens: `SIMPLE = T`, a size that is a multiple of 2,880,
  sane `BITPIX` and `NAXISn`, and every HDU's header and data accounted for by the size. A truncated
  download is fetched again;
- datacubes are unpacked there (`download.unpack_datacubes`), and the Astra transform is applied there, only to
  frames whose `OBSERVER` is Astra and that lack `ASTRAROT`/`ASTRAMIR` (ANDOR rotated, SPIRIT mirrored, by
  `download.astra_transform` as in the nightly download); every frame is verified again, and an Astra frame
  still without the keyword is not placed;
- each frame is hard-linked into the night directory only if no file with its stem exists there (any FITS
  suffix) and no differently named frame has its `DATE-OBS`; nothing is overwritten, the night directory is
  created only when there is a frame to put in it, nothing outside staging is deleted, and `download_log.csv`
  is not touched;
- every file added goes into `<queue_root>/fetch/manifest.csv`: time, telescope, night, action, file, ESO dp_id
  (the cube's, for unpacked frames), dp_type, object, bytes, md5 and job, so the additions can be audited or
  undone.

Two modes swap files instead of only adding, and only when asked: `--replace invalid` (`q add --repair`) also
re-fetches frames on disk that fail verification, such as the ~300 truncated raw frames (decision G3);
`--replace all` (`q add --refetch`) fetches the whole night again and swaps every frame that differs from the
fresh copy. A swap first hard-links the old file into `<queue_root>/quarantine/<TEL>/<NIGHT>/`, then renames
the verified fresh copy over it, so it is atomic and the old file is kept (the manifest records `replaced`,
the old size and md5, and where it went). A `.fits.fz` frame is never swapped for a `.fits` one.

A fetch job with `--rerun-class P1` (the look-back) queues a rerun of the night if it added science frames;
with `--rerun-class none` it only downloads, and whoever asked for it queues the pipeline run (below). Exit
0 means the pass finished (rows ESO does not have, or keeps serving broken, are listed in the log and the
result file); 75 means some rows failed for network reasons and may work later.

To try it without touching the archive, point a scratch queue root at a scratch archive: a `config.json` there
with `{"basedir": "<scratch>/base"}`, then `ORCHARD_QUEUE_ROOT=<scratch>/queue python -m jobqueue jobs fetch
--telescope Callisto --night 20261005 --rows <rows.json> --rerun-class none`, run from a directory holding the
code to test (`python -m` finds `jobqueue` in the working directory before `PYTHONPATH`). `--dest <dir>`
instead downloads into a plain directory and records nothing.

## Downloading and processing old nights

`q add --download` asks for both at once: an add-only fetch of whatever ESO holds for the night that the disk
lacks, then the pipeline run, as a job that waits for the fetch and is cancelled if the fetch fails. The ESO
query and the comparison run when the request is made, so `q add` says how many frames, and roughly how many
GB, the fetch will bring, and queues no fetch if nothing is missing. `--download-only` queues just the fetch.

A file of requests in the format of `manual_runs/backlog/v3_backlog.txt` (`TEL,DATE,DOWNLOAD,DELETE,PROCESS[,TARGETS]`)
queues them all: `q add --class P2 --from-file backlog.csv`. Every line is checked first, and nothing is queued
if any line is bad (an unknown telescope, a date that is not `YYYYMMDD`, flags other than 0/1, DELETE without
DOWNLOAD, DOWNLOAD for Artemis). `DELETE=1` is taken as `--refetch`: nothing is deleted. Compared with
`v3_manual_process.sh`, the same lines get the queue's priorities, locks, disk throttling, retries, per-job
logs and exit codes, and a failed download no longer leads to a pipeline run on a half-empty night.

The ESO backfill (decision D3) is the same thing with `--download-only` lines: what it adds is in the
manifest, and the nights it fills can then be queued for processing with `PROCESS=1` lines.

## The command line

Inside the container it is `python -m jobqueue <command>`; on appct, `src/jobqueue/cron/q` does the
`docker exec` for you (`q --shadow …` for the shadow database).

```bash
q status                       # dispatcher, classes, running and queued jobs (with what each waits for), failures
q status --resources           # also load, memory and disk %util (samples for 3 s)
q show 42                      # one job and all its events
q add --class P2 --telescope Europa --night 20250101          # a pipeline night (targets read from the frames)
q add --class P2 --telescope Europa --night 20250101 --no-T12 --targets "Sp0246+1625"
q add --class P2 --telescope Europa --night 20250412 --download         # fetch what is missing, then run it
q add --class P2 --telescope Europa --night 20251126 --download-only --repair   # also re-fetch truncated frames
q add --class P2 --from-file backlog.csv --dry-run                         # TEL,DATE,DOWNLOAD,DELETE,PROCESS[,TARGETS]
q add --class P3 --kind sweep --cores 2 --no-disk-heavy --estimate 30 -- python tools/x.py
q pause P2                     # resume P2 | pause all | resume all
q drain                        # start nothing new, let running jobs finish (q undrain)
q cancel 42                    # a queued job; a running one only with --kill (SIGTERM, then SIGKILL)
q signatures                   # failed jobs grouped by failure signature
q requeue --signature "gaia_dr3_id" --dry-run                 # then without --dry-run
q stop                         # stop the dispatcher; jobs carry on, the watchdog restarts it
q probe                        # read-only: load, memory, disk %util, code version, locks held outside the queue
q lookback --dry-run           # what the look-back would do today, queueing nothing
q report --days 3              # the daily report (live: jobs, results, delays; shadow: next to the cron)
```

## Shadow mode

The dispatcher, watcher and look-back run exactly as they would live, against their own database under
`mh_scratch`, while today's cron keeps doing the work. The shadow dispatcher launches nothing: an admitted job is
recorded as started and finishes after its estimate, so shadow jobs hold cores and locks realistically and the
admission decisions use appct's real load, memory and disk. Nothing in the production tree changes.

Every morning `queue-cron.sh --shadow report` writes `shadow/reports/shadow_<date>.md`: per telescope-night, when
the watcher found the night ready, when the pipeline would have started and finished, the cron's actual download
(start, end, attempts), its counts from `download_log.csv`, and its actual pipeline run (start and end from the
`date` lines of `<NIGHT>_1_v3.log`, in v3 or v2), with the gain in hours. Then jobs that had to wait and why,
the look-back's actions, and the watcher's notes.

### Starting it

A pinned snapshot of the queue code goes to `mh_scratch`, not the production tree. From a checkout of this
repository:

```bash
git archive <commit> src/jobqueue | ssh speculoos@appcs.ra.phy.cam.ac.uk \
    'mkdir -p /appct/data/SPECULOOSPipeline/mh_scratch/appct-job-queue && tar -x -C /appct/data/SPECULOOSPipeline/mh_scratch/appct-job-queue'
```

Then on appct as `speculoos`:

```bash
H=/export/data/SPECULOOSPipeline/mh_scratch/appct-job-queue
mkdir -p $H/shadow
echo '{"mode": "shadow", "queue_root": "/data/SPECULOOSPipeline/mh_scratch/appct-job-queue/shadow"}' > $H/shadow/config.json
crontab -l > ~/crontab.backup.$(date +%Y%m%d)
(cat ~/crontab.backup.$(date +%Y%m%d); grep -v '^#' $H/src/jobqueue/cron/crontab.shadow | grep .) | crontab -
```

Check it after 5 minutes with `q --shadow status`, and read the logs in `shadow/logs/`.

### Stopping it

```bash
crontab ~/crontab.backup.<date>       # the crontab as it was
q --shadow stop                       # the shadow dispatcher (it launches nothing, so there is nothing else)
```

## Going live (Q3, Q4)

After a clean week of shadow reports. On appct as `speculoos`, in the daytime:

1. **Deploy the queue code** with the rsync recipe from 7 October (stage `git archive <tag> src/jobqueue` on
   appct, then `rsync -rl --checksum --itemize-changes --backup --backup-dir=<dir>` into
   `/appct/data/speculoos/orchard/src/`, no `--delete`, no `-t`/`-p`). Only `src/jobqueue/`: it is a new
   directory, so nothing the pipeline runs changes. Check with `docker exec orchard-server python -m jobqueue
   probe`. Writing `src/VERSION` (decision K6) belongs with the next full deploy; until then jobs record a
   `tree-<hash>` fingerprint.
2. **Crontab**: save it (`crontab -l > ~/crontab.backup.<date>`), add the four queue lines from
   `src/jobqueue/cron/crontab.live` (watchdog, watcher, look-back, daily report), remove any shadow lines, and
   comment out the five nightly lines (19:00 Io … 23:00 Artemis). Leave the 05:00 and 06:00 lines.
   `crontab.live` is the whole proposed crontab, for reference.
3. **Watch the first day**: `q status`, `queue/logs/watcher.log`, the P0 jobs' logs, `ESO_logs/`, and the usual
   pipeline emails, which now arrive in the afternoon. Io's nightly INCOMPLETE email stops. Every morning
   `queue/reports/report_<date>.md` lists, per telescope-night, when it was ready at ESO, when its download and
   pipeline ran, how they ended, whether light curves reached v2 and the delay from ready to done, then failed
   jobs, look-back jobs and watcher notes. The queue sends no email of its own, so a broken watcher or dispatcher
   shows up only there and in `q status`.
4. The first live look-back may queue a dozen P1 jobs (see the dry run below). Run `q lookback --dry-run` first if
   you want to see them, or `q pause P1` until you have.

**Rollback:** `crontab ~/crontab.backup.<date>` restores the five nightly lines and removes the queue lines;
`q drain` stops new starts, and running jobs finish (or `q cancel --kill`); `q stop` stops the dispatcher. The
deployed `src/jobqueue/` directory can stay: nothing imports it.

## Configuration

Defaults are in `src/jobqueue/config.py`; a JSON file (`--config`, `$ORCHARD_QUEUE_CONFIG`, or
`<queue_root>/config.json`) overrides any of them, nested keys merged. `$ORCHARD_QUEUE_ROOT` overrides
`queue_root`. The ones most likely to change:

| Key | Default | |
|---|---|---|
| `mode` | `live` | `shadow` records and simulates only |
| `queue_root` | `/data/SPECULOOSPipeline/queue` | database, logs, staging, fetch lists, reports |
| `max_cores`, `load_max`, `mem_min_gb`, `disk_util_max` | 112, 100, 64, 70 | admission |
| `p0_reserve_cores`, `p0_reserve_window_utc` | 40, 11:00–01:00 | |
| `telescope_lock` | true | false once the per-run catcache lands |
| `semaphores.eso`, `semaphores.eso_bulk` | 2, 1 | concurrent download/fetch jobs; P2/P3 fetches at once |
| `fetch.workers`, `fetch.attempts`, `fetch.backoff_seconds`, `fetch.breaker` | 3, 4, 30, 8 | downloads in flight per job, tries per file, first backoff, network failures before a job gives up for now |
| `fetch.batch`, `fetch.mb_per_s` | 25, 4 | rows per staging batch; download rate assumed for estimates |
| `timeout_factor`, `min_timeout_minutes` | 3, 120 | |
| `retry_delay_minutes`, `max_retries` | 60, 2 | transient failures |
| `pipeline_cores` | 20 | sets `N_CORES` for the job |
| `minutes_per_1000_frames` | ANDOR 69, SPIRIT 22 | estimates (2026 production medians) |
| `watcher.*` | see file | stability, deadlines, `recent_days` = 5 |
| `lookback.*` | 30 nights, 2 reruns/night | |
| `p0_guard.windows_utc` | SSO 13:00–17:00, Artemis 10:00–13:00 | |

## Files

| | |
|---|---|
| `src/jobqueue/store.py` | SQLite store and state transitions |
| `src/jobqueue/scheduler.py` | admission and strict priority (pure) |
| `src/jobqueue/dispatcher.py` | the daemon, watchdog, code version |
| `src/jobqueue/runner.py` | runs one job, writes its status |
| `src/jobqueue/resources.py` | `/proc`: load, memory, disk %util, liveness, locks held outside the queue |
| `src/jobqueue/signatures.py` | failure signatures, transient or deterministic |
| `src/jobqueue/eso.py` | authenticated read-only TAP |
| `src/jobqueue/nights.py` | disk side: cube-aware diff, FITS headers, products, transfer log |
| `src/jobqueue/watcher.py`, `lookback.py`, `jobs.py`, `jobspec.py` | nightly, look-back, job bodies and descriptions |
| `src/jobqueue/fetch.py` | the add-only fetch: verified, paced downloads, placing, swaps, manifest |
| `src/jobqueue/shadow_report.py` | the shadow report |
| `src/jobqueue/cli.py` | `q` |
| `src/jobqueue/cron/` | `queue-cron.sh`, `q`, `crontab.shadow`, `crontab.live` |
| `src/jobqueue/tests/` | `pytest src/jobqueue/tests` |

## Limits and next steps

- **One pipeline job per telescope** until the per-run catcache (K7). Then set `telescope_lock: false`.
- **Emails.** P0 jobs send what the cron sends today. P1 reruns send the normal pipeline emails until a quiet
  switch exists (K7).
- **No product check for P2 yet.** The plan's verify-then-promote step (run with `--no_T12`, check the ledger,
  then promote and touch) belongs with the P2 backlog work.
- **`SSO_download.py` with no transfer-log entry** reports "INCOMPLETE" in its email even when the night is
  complete, because its success test needs the transfer log.
- **Other containers** are invisible to the external-lock scan.
- Queue depth, throughput and failures in the nightly ledger report (R1) are not wired up yet; `q status` and
  `q signatures` cover it meanwhile.
