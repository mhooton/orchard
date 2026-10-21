# The appct job queue

A priority queue for the orchard pipeline on appct (`src/jobqueue/`), an ESO watcher that queues each night as
soon as its frames are at ESO, and a daily look-back over the last 30 nights. It is decision Q1 of the portal
completeness plan, with Q3 and Q4 (nightly processing through the queue, download separate from processing) as
the first use. Backlog reruns (class P2) come later, once the 16-bit frames, a per-run catcache (done, K7) and the
K8 fixes have landed.

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

- declared cores of running jobs plus its own fit in **112** of the 128 threads (`p0_reserve_cores` can keep
  some of them for P0 from 11:00 to 01:00 UTC; it is 0, so nothing is held back);
- the 5-minute load average, plus the cores started earlier in the same poll, is under **100**;
- `MemAvailable` is over **64 GB**;
- both arrays (`sda` = PipelineOutput, `sdb` = data/raw) are under **70 %** utilisation, measured from
  `/proc/diskstats` over the last poll interval. P0 skips this check. At most one disk-heavy job starts per
  poll, because `%util` only shows a new job's I/O at the next poll;
- its locks are free (below).

Several nights of one telescope may run at once. Nothing is kept free for a nightly job that has not arrived:
when last night's data is ready it goes to the front of the queue and takes the next capacity that frees up,
which may mean waiting for a running backlog job to finish. (Until 2026-10-17 a `p0_guard` kept backlog jobs
off a telescope from 13:00 to 17:00 UTC, and 40 cores were kept for P0.)

### Strict priority, never preemption

Jobs are considered by class, then priority within the class, then age. When a job cannot start, the resources
it is waiting for are reserved for its class: no lower-class job may take them. A lower-class job may still use
something else (a P1 job runs while a P0 job waits for a target another job holds), and a same-class job may use
what is left (a P0 download is not held up by a P0 pipeline job waiting for cores).

A running job is never killed or paused to make room: the plate solver relies on a 60 s `SIGALRM`, and the
pipeline holds files open. A P0 job that finds the machine full waits for the next job to finish. The only kills are a job
past three times its estimate (floor 120 min) and an explicit `q cancel --kill`.

### Locks

| Lock | Held by | Why |
|---|---|---|
| `tel:<TEL>` | pipeline jobs | Enforced only with `telescope_lock: true` (one pipeline job per telescope). Off since K7, which removed what two nights of one telescope shared: each wcsfit call has its own catcache; `reports/Pipeline_Report`, `reports/QC`, `output/StackImages/already_run.txt` and each star's `output/<gaia_id>_<filter>_output.fits` are updated under a file lock (`utils/shared_files.py`); the PDF report (T11) works in `reports/temp/<night>`. |
| `night:<TEL>:<NIGHT>` | every job on that night | raw and output directories of one night |
| `target:<NAME>` | pipeline jobs, per target | `PipelineOutput/v2/StackImages` is shared across telescopes and nights |
| `sem:eso` | download and fetch jobs | at most `semaphores.eso` (2) at once |
| `sem:eso_bulk` | P2 and P3 fetch jobs as well | at most `semaphores.eso_bulk` (1), so bulk downloads never hold both ESO slots and the nightly download always finds one free |

Pipeline and download processes the queue did not start (the cron, a person, a replay) hold the same locks: the
dispatcher reads them from `/proc/*/cmdline` inside the container. So does a promotion run that changes the tree
(`utils/promote.py ... --apply --telescope T --night N`), which holds that night. Processes in other containers
(`orchard-server-spirit`, `orchard-server-autosave`) are invisible to it.

Targets for locking come from ESO's `object` column (SSO) or the raw file names (Artemis), upper-cased with
spaces as `--`, as `createlists.py` names them.

### Jobs, logs, retries

Each job runs through `jobqueue.runner` in its own session, so the job and everything it starts share one
process group. Output goes to `<queue_root>/logs/jobs/<YYYY-MM>/<id>-<kind>-<tel>-<night>-a<attempt>.log`, which
unlike `PipelineOutput/v3/<TEL>/logs/<NIGHT>_1_v3.log` is never overwritten by a later run. The runner writes
`<log>.status.json` with the exit code when the command ends.

`<attempt>` counts every run the job has had, a requeue included, and the dispatcher never gives a run a log or
status file that already exists, so each run has files of its own. A status file counts only if the runner
that wrote it started no earlier than the attempt (60 s of slack for a clock step); an older one is ignored,
with a `stale-status` event. Before this, `q requeue` reset the numbering with the retry count. On 19 October
jobs 1502, 1625, 1627 and 1649 got attempt 1's paths again. The next poll failed each one on attempt 1's exit
status and signature, freeing its locks and cores, while the run carried on and ended with exit 0.

Every job records its start, end, exit code and the code it ran: the release it was pinned to (see *Deploying
code*), or, with no release, the first word of `src/VERSION` if it exists, otherwise `tree-<hash>` of the
deployed `.py` and `.sh` files (the production tree is a working copy,
not a git checkout).

When a job fails, its **failure signature** is the last exception line of its log (else the last error line),
with paths, numbers and target names normalised. It is **transient** if the signature or the last 30 lines
match a network, full-disk, I/O or memory pattern, or the exit code is 75, 137 or 143; transient failures are
retried after an hour, at most twice. Anything else is **deterministic**: the job stops as failed, and after a
fix `q requeue --signature "<text>"` puts every job with that signature back, with its retries restored.

### Threads

orchard-server's systemd scope (`docker-<id>.scope`, systemd 219 on CentOS 7) allows **4,096 tasks for the whole
container**, and every thread counts. Docker sets no limit of its own (`PidsLimit` 0). By default every Python process
that uses numpy starts OpenBLAS's pool of up to 64 threads, so one T7 with 20 workers held 1,280.

So jobs run with **one BLAS/OpenMP thread per process**. `jobqueue.runner` sets `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS` and `MKL_NUM_THREADS` to 1 for every job, unless the job's environment sets them, and
`ZLP_pipeline.sh` does the same for runs outside the queue. Its log says `BLAS threads per process: ...`. `q release
test` gives the suite the same. Orchard's parallelism is its process pools, and its BLAS calls are 2×2 inverses and
small least-squares fits, so a single thread costs nothing.

**What exhaustion looks like** (2026-10-20/21, v4.1.0 to v4.3.2: the new T7 solver's `np.linalg.inv` started the pool
in every worker). When `pthread_create` fails, OpenBLAS prints `OpenBLAS blas_thread_init: pthread_create failed for
thread N of 64` and raises SIGINT:

- At numpy import, the stage dies with `KeyboardInterrupt`, exit 130, within a minute of the start. The queue counts it
  as transient, and the retries are spent.
- In a pool worker, the worker dies mid-chunk (`Process ForkPoolWorker-N: ... KeyboardInterrupt`). `Pool.map` waits for
  the lost chunk for ever, and the job sits until its timeout. Its idle workers keep their 64 threads, which makes the
  shortage worse.
- Other signatures from the same cause: `RuntimeError: can't start new thread`; `AttributeError: 'DummyProcess' object
  has no attribute 'terminate'` (a `ThreadPool` failing to start; `ThreadPool()` in createlists, pipered and condense
  starts one thread per core, i.e. 128); `ImportError: numpy._core.multiarray failed to import`.

To check, on appct run `systemctl show docker-$(docker inspect --format '{{.Id}}' orchard-server).scope -p TasksMax
-p TasksCurrent`, or inside the container read `/sys/fs/cgroup/pids/pids.current` and `pids.max`. Scratch work in the
container shares the same 4,096: give it the same three variables.

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
(pauses, drain, dispatcher heartbeat, `eso_reach`: whether ESO answered the latest polls), and the notifier's
`notify_incidents`, `notify_messages` and `job_pwv` (the PWV outcome of each finished SSO pipeline job).

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
- If ESO cannot be reached the watcher changes nothing; the next poll retries. Each poll (and each look-back run)
  records whether ESO answered any of its queries in `kv` `eso_reach`, and the notifier opens an incident once
  the failures have lasted 4 polls or 2 hours of watcher time (below).

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
  their copies take twice the disk). A row whose file sits in the previous or next night's directory is
  present too: 161 frames, almost all morning darks, were filed under the next night;
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
DOWNLOAD, DOWNLOAD for Artemis). `DELETE=1` is taken as `--refetch`: nothing is deleted. A PROCESS line without TARGETS
has its targets (for the locks) read from the night's frames when it is queued, which for an SSO night means the
header of every raw frame. `--dry-run` skips that, but queueing a long file does not: about 1.5 million frame
headers for a 1,156-night file. Give TARGETS for long files. Compared with
`v3_manual_process.sh`, the same lines get the queue's priorities, locks, disk throttling, retries, per-job
logs and exit codes, and a failed download no longer leads to a pipeline run on a half-empty night.

The ESO backfill (decision D3) is the same thing with `--download-only` lines: what it adds is in the
manifest, and the nights it fills can then be queued for processing with `PROCESS=1` lines.

## Promotion after pipeline jobs

T12 is being replaced by the promotion step, `utils/promote.py` ([promotion.md](promotion.md)). A pipeline job queued
with `--promote` runs with `--no_T12`, and a **promote** job follows it:

- It runs `python utils/promote.py run --apply --telescope T --night N --job <pipeline job>`. That judges the targets
  the pipeline job ran, using `promote_check`. What passes is published into v2, and the rest is kept in the
  promotion journal (`<queue_root>/promote.sqlite`) for a person.
- It runs whether the pipeline job succeeded or failed, so the targets that finished are judged. A cancelled
  pipeline job cancels it.
- It holds `night:<TEL>:<NIGHT>`, takes 2 cores, is not disk-heavy, and runs ahead of the pipeline jobs of its class.
- `q requeue` of the pipeline job requeues it as well.

The daily report lists the holds and rejects waiting for a decision, and any promotion that did not finish. The
nightly runs from the watcher and the look-back still use T12 until `promote.nightly` is turned on.

**Reruns from a late stage** (`q add --from T10 --promote`, or `--from pwv`; [rerun-from-stage.md](rerun-from-stage.md))
queue a **restage** job first (`utils/restage.py`: hard links of the rerun's inputs from v2 into v3, 1 core, the night
lock, a second or two), then the pipeline job with `--from_T10 --no_T12` for the restaged targets, which runs only if
the restage succeeded, then the promote job, which overlays v2 with what the rerun wrote.

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
q add --class P2 --from-file backlog.csv --promote --campaign "ESO backfill reruns"   # --no_T12, then a promote job each
q add --class P2 --telescope Europa --night 20251016 --promote             # one night the same way
q add --class P2 --kind promote --after 1544-1670                          # promote jobs after pipeline jobs already queued
q add --class P2 --kind promote --telescope Europa --night 20251016        # judge one night against its latest run
q add --class P3 --telescope Europa --night 20260703 --from pwv --promote [--dry-run]   # restage, rerun T10+T11, overlay
q add --class P3 --from-file nights.csv --from T10 --promote --campaign "PWV outage reruns"  # PROCESS-only lines
q promote list                 # holds and rejects waiting for a decision (q promote = python utils/promote.py)
q promote show 12              # q promote accept 12 --reason "..." [--apply]; q promote reject 12 --reason "..."
q add --class P3 --kind sweep --cores 2 --no-disk-heavy --estimate 30 -- python tools/x.py
q pause P2                     # resume P2 | pause all | resume all
q drain                        # start nothing new, let running jobs finish (q undrain)
q cancel 42                    # a queued job; a running one only with --kill (SIGTERM, then SIGKILL)
q signatures                   # failed jobs grouped by failure signature
q requeue --signature "gaia_dr3_id" --dry-run                 # then without --dry-run
q stop                         # stop the dispatcher; jobs carry on, the watchdog restarts it
q probe                        # read-only: load, memory, disk %util, code version, locks held outside the queue
q release                      # code releases; q release deploy --from <stage>/src --sha <sha> (see Deploying code)
q lookback --dry-run           # what the look-back would do today, queueing nothing
q report --days 3              # the daily report (live: jobs, results, delays; shadow: next to the cron)
q add --class P3 --telescope Callisto --night 20250301 --campaign "WCS audit reruns"   # grouped so in the digest
q campaign "WCS audit reruns" --note "WCS audit" --dry-run    # label jobs already queued (or --ids 1200-1350)
q notify digest --print        # what this morning's digest would say; changes and delivers nothing
q notify incidents [--all]     # open incidents (and recently resolved ones)
q notify ack 7 --until 2026-11-01 --note "dome repair booked"   # hold its alert; it stays in the digest
q notify messages              # what was delivered, with the outbox file of each
q notify recipients --add cj467@cam.ac.uk     # the recipient list (notify.to in queue/config.json); --remove
```

## Notifications

`notify.py` is the pipeline's only sender of email; `digest.py` decides what it says. Pipeline stages should not
email anyone; `ZLP_pipeline.sh` stopped calling `email_eon.py` once the digest was sending (18 Oct 2026), and the
script was deleted on 19 Oct 2026.

- **Incidents.** Every run (`queue-cron.sh alerts` every 30 minutes, `queue-cron.sh digest` each morning) works out
  the conditions that need a person and keeps one incident per condition in the queue database
  (`notify_incidents`). An incident gets one email when it opens and one when it resolves, in the same thread;
  in between it is listed in each digest. `q notify ack` holds the alert until a date. The conditions:

  | Key | Opens when | Resolves when |
  |---|---|---|
  | `failed:<tel>:<night>:<kind>` | the newest download or pipeline job of a night that had a P0 job failed (last 3 nights) | a later job of the same kind succeeds |
  | `late:<tel>:<night>` | the transfer log lists frames ESO still lacks at 21:00 UTC the next day, and the watcher is still waiting | the watcher queues the night |
  | `dark:<tel>` | 3 nights with no survey data on which two other SSO telescopes observed (any 3 nights for Artemis) | the telescope observes again |
  | `signature:<hash>` | 10 failures with one signature in 24 h, any class | fewer than 10 in the last 24 h |
  | `dispatcher`, `queue-stalled` | no dispatcher heartbeat for 30 min; jobs due but none running or started for 6 h | heartbeat back; a job starts |
  | `disk:<label>` | a disk under 2 TiB or 5% free | space is back |
  | `eso-unreachable` | every ESO query of the watcher's last 4 polls failed, or failed for 2 h of watcher time (look-back runs count as polls; the overnight break does not count as time) | ESO answers a query |
  | `pwv-nodata` | 3 SSO pipeline jobs since the last one that corrected PWV logged `No PWV data for:` on most of their light curves | a later job corrects PWV |

  The last two come from the outage of 18 October 2026, when appct's first nameserver began answering SERVFAIL,
  nothing resolved inside the containers, and for a day and a half the watcher logged a failed ESO query every
  half hour while 95 pipeline jobs wrote light curves without PWV. `eso_reach` holds the run of failed polls (from
  `watcher.py` and `lookback.py` through `esoreach.py`; a watcher that cannot set up its ESO client counts too).
  For PWV, `pwvlog.py` reads each finished SSO pipeline job's log once, newest first and at most 40 a run, for jobs
  finished in the last 48 h, and keeps the counts in `job_pwv`: each `python …SPlightcurve.py` line is one light
  curve, and `PWV correction successful!`, `No PWV data for:` or `PWV Failed: No T_eff for target` (off by design:
  Artemis, Teff < 2000 K) says what became of it. The incident names the jobs and the cause the log gives, and it
  stays open while no SSO job runs, until one corrects PWV.

  A night with no frames is attributed, not explained: a night SPOCK planned mostly as Chilean time (`ch_` blocks
  in `Archive_night_blocks`) has no survey frames by design and never counts; a telescope dark while two others
  observe is down itself; when none observes it is the weather, which needs nobody.
- **The digest** of last night (at 07:45 the night before last, processed yesterday) goes out every morning even
  when all is well: one line per telescope (frames, hours, light curves, or what went wrong, and the frames the
  telescope itself lost as the night report recorded them: `docs/telescope-faults.md`), open and newly
  resolved incidents, messages that could not be sent, failure signatures, look-back actions, ESO query errors,
  the PWV outcome of the SSO light curves made in the last 24 h (corrected, without LHATPRO data, other failures,
  off by design), free disk, and the queue's work in two sections.
  *Processed in the last 24 h* and *Remaining* group jobs by class and campaign: done, failed and cancelled
  counts, telescopes, the nights covered, light curves made (for pipeline work), estimated hours left, and an
  ETA at the last 24 hours' pace. A job's campaign is its `q add --campaign` label (`q campaign` labels jobs
  already queued); without one, P0 is "nightly", P1 "look-back", and anything else is grouped by its note.
  If it stops arriving, something upstream of it is broken; `notify.heartbeat_url` (fetched after each digest)
  lets an outside service such as healthchecks.io say so when the whole server is down.
- **Delivery.** Every message is written to `queue/notify/outbox/<time>_<kind>-<subject>.eml` and recorded in
  `notify_messages`. Nothing is sent until `notify.transport` is `smtp` in `queue/config.json`; then each message
  also goes to the recipient list, `notify.to` (`q notify recipients` shows it, `--add` / `--remove` edit it in
  place), with `EMAIL_ADDRESS`, `EMAIL_PASSWORD` and `SMTP_SERVER` read from `notify.env_file`, a file readable
  only by speculoos (the live `reporting/.env`, so a new Gmail app password reaches the notifier too). A message
  that cannot be sent still lands in the outbox; an alert that failed is tried again on the next run, and the next
  digest that gets through lists every message since the last delivered digest whose sending failed, with the
  error. The SMTP host is a name, so while appct cannot resolve names (as on 18 October 2026) nothing is sent, and
  the outbox and that later digest are the record. Incidents already alerted to the outbox before the switch are
  not re-sent; the digest lists them.

  ```json
  {"notify": {"transport": "smtp", "env_file": "/opt/orchard/src/reporting/.env",
              "to": ["mh2143@cam.ac.uk"], "heartbeat_url": null}}
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
   `crontab.live` is now a copy of the whole crontab installed on appct, for reference.
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

## Deploying code: releases

Every job runs from a **release**: a frozen copy of the production tree, so a deploy never waits for a gap in
the queue. Before releases, a deploy had to wait until no `ZLP_pipeline.sh` was running, because a night that
ran T8 on the old code and T9 on the new would mix them; with three telescopes and a backlog that gap could be
days away, and on 16 October a Callisto night was killed to make one.

- A release is `<queue_root>/releases/<id>/src` (on appct `/export/data/SPECULOOSPipeline/queue/releases`), with
  `<id>` = `<UTC time>-<sha7>`. It is the live tree (`/opt/orchard/src`) as it was when the release was built,
  with the staged commit laid over it, so deployed-only files (the target tables) come along.
  Caches, logs and core dumps in the live tree stay out (`release_exclude`). Its files are read-only and keep
  their own read bits. Files unchanged since the previous release, mode included, are hard links to it, so a
  release costs only what changed.
- Credentials are not copied (`release_link`, by default every `.env`): for `reporting/.env` a release holds a
  symlink to `/opt/orchard/src/reporting/.env` instead. The file keeps its own mode and exists only once, and a
  rotated password reaches every release at once. Neither a live-tree sync nor a rollback writes one: releases built before this change hold a
  copy, and activating one of those leaves the live file alone (its jobs still read its own copy). The symlink
  resolves only inside the container.
- `releases/current` is a relative symlink to the release new jobs get. When the dispatcher starts a job it
  points the job's working directory, `ORCHARD_PATH` and the live-tree entries of `PYTHONPATH` at that release,
  so `ZLP_pipeline.sh` (`SCRIPTDIR=$ORCHARD_PATH`), every script it runs and `python -m jobqueue jobs …` read the
  release from the first stage to the last. Activating another release changes nothing for running jobs; a retry
  starts on whatever is current then. Paths outside the live tree (a scratch checkout) are left alone.
- The job's `code_sha` is the release id, `q status` shows the current release and marks running jobs on an
  older one, and each job log's header has its `ORCHARD_PATH`.
- The live tree is still what everything outside the queue runs: the dispatcher, the watcher and look-back,
  cron, `q`, and people. Activating a release marks the live tree for an update, which writes the files that
  differ (each through a temporary file and a rename, as the rsync recipe did), backs up what it replaces in
  `releases/live_backups/`, and deletes nothing. It waits while a pipeline, download or job process still runs
  from the live tree (one started outside the queue, or by the queue before releases existed), and the
  dispatcher does it on the first poll after the last one ends. Pinned jobs and anyone's scratch copy (an
  `ORCHARD_PATH` elsewhere) do not hold it up. If the update rewrites the queue's own code, the dispatcher exits
  and the watchdog starts the new one within 5 minutes; running jobs carry on.
- With no `current` release the dispatcher launches from the live tree exactly as before.

### Deploying main

If the deploy changes products, tag the commit first (`git tag -a v4.1.0 -m ... <sha>; git push origin v4.1.0`)
so its products carry a version number of their own; see [versioning.md](versioning.md). Stage the commit with
`git archive`, which writes the version into the staged tree, then one command on appct does the rest:

```bash
git archive <sha> src | ssh speculoos@appcs.ra.phy.cam.ac.uk \
    'mkdir -p /appct/data/speculoos/deploy_staging/main_<sha> && tar -x -C /appct/data/speculoos/deploy_staging/main_<sha>'
q release deploy --from /gaia_database/deploy_staging/main_<sha>/src --sha <sha>
```

`deploy` takes a minute or two whatever the queue is doing. It builds the release, runs the suite from it
(`tests/` and `jobqueue/tests/`, with the pytest pylib in `test_pylib`; it stops before activating if anything
fails, and `--skip-tests` skips it), activates it, and brings the live tree up to it straight away if nothing
runs from the live tree, else says what it is waiting for and leaves the update to the dispatcher. If the
queue's own code changed it restarts the dispatcher (`q stop`; the watchdog starts the new one within 5
minutes; running jobs carry on). The suite's temporary files (pytest's `--basetemp` and `TMPDIR`) go to
`<release>/pytest-tmp` on the data disk and are removed when it ends, pass or fail; `tests.log` and `tests.json`
stay in the release. The steps also run one at a time:

```bash
q release build --from <stage>/src --sha <sha>   # or no --from: a snapshot of the live tree as it is
q release test <id>
q release activate <id>                          # jobs started from now on run it; the live tree follows
q release                                        # releases, which jobs run on each, whether the live tree matches
q release sync-live --dry-run                    # what differs in the live tree
q release sync-live [--force]                    # update the live tree now (--force: even while it is in use)
```

**Rollback:** `q release activate <previous id>`. New jobs get the old code at once; the live tree follows.

**When to hold the queue anyway:** a change that must not overlap with runs of the old code (a new output
layout two runs would both write, say) still needs `q drain`, waiting for the running jobs, then the deploy and
`q undrain`. That is the exception; a fix to a pipeline stage is not.

**Data files in the tree** (`calibration/*.csv`, `utils/ml_40pc_v2.csv`) reach jobs with the next release. If
one is edited in the live tree, snapshot it with `q release deploy` (no `--from`). The same goes for any
hand edit in the live tree: until the next release, jobs do not see it.

### Running the suite in orchard-server yourself

**Always give pytest `--basetemp` and `TMPDIR` on `/data`, and delete them afterwards.** The container's root,
`/tmp` included, is the overlay on appct's 20 GB `/var`, and `df` inside the container shows `/var`, not the data
disk. A full suite leaves about 1.1 GB in `/tmp/pytest-of-unknown` and pytest keeps the last three runs. On
2026-10-19 such runs filled `/var` twice, and production jobs 1877 and 1891 died with `Errno 28` in T9 and T10.
The container has no pytest of its own: use the pylib, and run from the tree under test (`python -m` puts the
working directory first on `sys.path`, and the container's default is the live tree). On appct, for a
`git archive` stage in `mh_scratch/<dir>`:

```bash
docker exec orchard-server bash -c '
  T=/data/SPECULOOSPipeline/mh_scratch/<dir>; P=$T/pytmp
  mkdir -p $P/tmp && cd $T/src &&
  TMPDIR=$P/tmp PYTHONPATH=/data/SPECULOOSPipeline/mh_scratch/k1-deploy-check/pylib:$T/src PYTHONDONTWRITEBYTECODE=1 \
    nice -n 15 python -m pytest -q -p no:cacheprovider --basetemp=$P/basetemp tests jobqueue/tests
  rm -rf $P'
```

`--basetemp` must be a directory of its own: pytest empties it when it starts. `q release test` and
`q release deploy` do all of this themselves.

## Configuration

Defaults are in `src/jobqueue/config.py`; a JSON file (`--config`, `$ORCHARD_QUEUE_CONFIG`, or
`<queue_root>/config.json`) overrides any of them, nested keys merged. `$ORCHARD_QUEUE_ROOT` overrides
`queue_root`. The ones most likely to change:

| Key | Default | |
|---|---|---|
| `mode` | `live` | `shadow` records and simulates only |
| `queue_root` | `/data/SPECULOOSPipeline/queue` | database, logs, staging, fetch lists, reports |
| `max_cores`, `load_max`, `mem_min_gb`, `disk_util_max` | 112, 100, 64, 70 | admission |
| `p0_reserve_cores`, `p0_reserve_window_utc` | 0, 11:00–01:00 | cores kept for P0 in that window |
| `telescope_lock` | false | true: one pipeline job per telescope, including runs started outside the queue |
| `semaphores.eso`, `semaphores.eso_bulk` | 2, 1 | concurrent download/fetch jobs; P2/P3 fetches at once |
| `fetch.workers`, `fetch.attempts`, `fetch.backoff_seconds`, `fetch.breaker` | 3, 4, 30, 8 | downloads in flight per job, tries per file, first backoff, network failures before a job gives up for now |
| `fetch.batch`, `fetch.mb_per_s` | 25, 4 | rows per staging batch; download rate assumed for estimates |
| `timeout_factor`, `min_timeout_minutes` | 3, 120 | |
| `retry_delay_minutes`, `max_retries` | 60, 2 | transient failures |
| `pipeline_cores` | 20 | sets `N_CORES` for the job |
| `releases_dir` | `<queue_root>/releases` | code releases; `current` is the one new jobs get |
| `release_exclude` | see file | caches, logs, core dumps left out of a release (`/x` anchors at the tree's root) |
| `test_pylib`, `test_paths` | `mh_scratch/k1-deploy-check/pylib`, `tests`, `jobqueue/tests` | `q release test` |
| `test_tmp` | `<release>/pytest-tmp` | the suite's `--basetemp` and `TMPDIR` (`<test_tmp>/<release id>` if set), removed after each run; never the container's `/tmp` |
| `minutes_per_1000_frames` | ANDOR 69, SPIRIT 22 | estimates (2026 production medians) |
| `telescopes.<TEL>.force_platesolve` | true for Io, Europa, Ganymede, Callisto; false for Artemis | pipeline jobs pass `--force-platesolve`: T7 ignores the night's WCS cache, re-solves every frame and refreshes the cache rows ([wcs-cache.md](wcs-cache.md)). SSO stays forced until the cache is validated |
| `watcher.*` | see file | stability, deadlines, `recent_days` = 5 |
| `lookback.*` | 30 nights, 2 reruns/night | |
| `notify.transport`, `notify.env_file`, `notify.to` | `outbox`, none, none | `smtp` also sends each message (see Notifications) |
| `notify.heartbeat_url` | none | fetched after each digest, for an outside dead-man's switch |
| `promote.db`, `promote.raw_cache` | `<queue_root>/promote.sqlite`, `<queue_root>/promote_raw.sqlite` | the promotion journal; raw-frame counts for check A4 |
| `promote.identity` | `plan` | A3 from plan files, or `table`: `<basedir>/ml_40pc_v2.csv` for Sp targets, holding where the plan disagrees |
| `promote.thresholds` | `{}` | overrides of `promote_check.THRESHOLDS`, e.g. `{"max_scatter_ratio": 1.3}` |
| `promote.cores`, `promote.est_minutes`, `promote.priority` | 2, 15, 1e9 | the promote job |
| `promote.nightly` | false | true: the watcher's, look-back's and fetch reruns' pipeline jobs run with `--no_T12` and a promote job |
| `notify.*` thresholds | see file | dark nights (3), late hour (21:00), stall (6 h), dispatcher (30 min), disk (2 TiB / 5%), signature (10/24 h), ESO (4 polls / 2 h), PWV (3 jobs; logs of the last 48 h, 40 a run) |

## Files

| | |
|---|---|
| `src/jobqueue/store.py` | SQLite store and state transitions |
| `src/jobqueue/scheduler.py` | admission and strict priority (pure) |
| `src/jobqueue/dispatcher.py` | the daemon, watchdog, code version |
| `src/jobqueue/runner.py` | runs one job, writes its status |
| `src/jobqueue/releases.py` | code releases: build, pin a job, activate, bring the live tree up to date |
| `src/jobqueue/resources.py` | `/proc`: load, memory, disk %util, liveness, locks held outside the queue |
| `src/jobqueue/signatures.py` | failure signatures, transient or deterministic |
| `src/jobqueue/eso.py` | authenticated read-only TAP |
| `src/jobqueue/nights.py` | disk side: cube-aware diff, FITS headers, products, transfer log |
| `src/jobqueue/watcher.py`, `lookback.py`, `jobs.py`, `jobspec.py` | nightly, look-back, job bodies and descriptions |
| `src/jobqueue/fetch.py` | the add-only fetch: verified, paced downloads, placing, swaps, manifest |
| `src/jobqueue/shadow_report.py` | the shadow report |
| `src/jobqueue/digest.py` | what the notifier says: per-night status, queue, disks, incident conditions, the digest text |
| `src/jobqueue/notify.py` | incidents, alerts and the digest; the outbox and SMTP delivery; the heartbeat |
| `src/jobqueue/esoreach.py` | whether ESO answered: each watcher poll's and look-back run's outcome, for the notifier |
| `src/jobqueue/pwvlog.py` | the PWV outcome of each finished SSO pipeline job, read from its log once |
| `src/jobqueue/cli.py` | `q` |
| `src/utils/promote.py`, `src/utils/promote_journal.py` | the promotion step and its journal ([promotion.md](promotion.md)) |
| `src/jobqueue/cron/` | `queue-cron.sh`, `q`, `crontab.shadow`, `crontab.live` |
| `src/jobqueue/tests/` | `pytest src/jobqueue/tests`; in orchard-server, see [Running the suite in orchard-server yourself](#running-the-suite-in-orchard-server-yourself) |

## Limits and next steps

- **Same-target nights stay one at a time**: the `target:` lock, because each target's `<gaia_id>_<filter>_output.fits`
  in `output/` gathers every night of that target on that telescope.
- **Emails.** The digest replaced the per-target `email_eon.py` emails, and the script was deleted on 19 Oct 2026, so
  `notify.py` is the only sender left in `src/`. The Pipeline_v2 senders are outside this repository.
- **Promotion is phase 1.** Pipeline jobs queued with `--promote` are judged and published by `utils/promote.py`;
  the nightly runs still use T12 (`promote.nightly`), and per-run staging and provenance keywords come after the
  backfill (see [promotion.md](promotion.md)).
- **`SSO_download.py` with no transfer-log entry** reports "INCOMPLETE" in its summary even when the night is
  complete, because its success test needs the transfer log.
- **Other containers** are invisible to the external-lock scan.
- Queue depth, throughput and failures now reach the morning digest; the nightly ledger itself (R1) is still to
  come, and so are notifier conditions that need it (images taken at the site but never archived, frames that
  vanish from disk).
