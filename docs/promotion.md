# Promotion: publishing v3 into v2 through the checker

`src/utils/promote.py` is phase 1 of the replacement for T12, the last stage of `ZLP_pipeline.sh`
(`perform_version_migration`). T12 moves pipeline output from `PipelineOutput/v3` into `PipelineOutput/v2`, the
tree the portal reads. The promotion step publishes from v3 into v2 only what the read-only checker
(`src/utils/promote_check.py`, PR #16) says is fit to replace what v2 has. Everything else waits in a hold queue for
a person. Pipeline jobs run with `--no_T12`, and a `promote` job of the queue runs this after each of them.

v2 stays where it is, because the portal hard-codes it.

## Why

T12 sweeps a v3 target directory into v2 whenever it holds a `*4_diff.fits`, and looks at nothing else. The
2026-10-09 audit (`mh_scratch/sweep_audit/`) found that it:

- nests on a night's second supersession (`mv` to a fixed `v2_old` path lands inside the first copy) and fails on
  the third (`File exists`, which kills the run inside T12);
- swaps the whole night's v2 reduction when only one target was rerun, so the other targets' products end up
  beside masters they were not made with;
- strands good results in v3 when a later target crashes, because it only runs at the end;
- made 13 products worse: 4 lost their aperture-5 curve, 2 now publish an all-NaN curve, 2 lost points and 5 got
  noisier.

The checker holds or rejects all 13 regressions. With R1 at 0.98 it holds 5 of the 422 good sweeps (about 1%), and
it promotes 86 of the 93 curves stranded in v3; the other 7 have no usable aperture 5.

## What is judged

A **candidate** is a v3 target directory, `v3/<TEL>/output/<NIGHT>/<TARGET>/`, that holds a `*_diff.fits`. Its
**incumbent** is the v2 directory at the same path, if there is one. The directory is replaced as a whole, so every
star in it is judged, not just the target.

The candidate is meant to be the output of a **run**: the latest pipeline run of the night that covers the target.
A run given targets (a single-target rerun) covers only those; a run without them covers the targets it
discovered. Runs come from two places:

- the queue's pipeline jobs (start and end of the latest attempt);
- the pipeline logs `<NIGHT>_1_v3[_<targets>].log` in `v3/<TEL>/logs` and `v2/<TEL>/logs`. T12 moves a night's log
  to v2 when it sweeps any target of it, and a later run of the night overwrites the log in v3.

A queue job is matched to its log where the log is still there. The run's start is its first `date` line, in UTC
inside the container.

Only a run that made light curves counts as the candidate's run: one that ran T10. A queue job run with
`--only_T11` (the PDF and frame packing, e.g. for a night whose run died after its science) or `--only_T12` makes
none, so S1 and the reduction are judged against the run before it. Such later runs are kept on the run as `later`,
and their PDF and log go into the records.

### Checks

From `promote_check.assess()`, unchanged except R1 (thresholds as agreed; every one can be overridden with
`--<name>`, or in the queue config under `promote.thresholds`):

| Check | Holds or rejects when | Threshold |
|---|---|---|
| A1 | no aperture-5 curve at all: REJECT; only for other stars: HOLD | the portal reads `*_5_diff.fits` only |
| A2 | no finite target point: REJECT; fewer than 20, or under half the frames: HOLD | `min_valid_points` 20, `min_valid_fraction` 0.5 |
| A3 | the curve's Gaia ID is not the plan's (aliases and the DR2/DR3 crosswalk applied): HOLD | |
| A4 | frames in the curve / raw science frames on disk below 0.8: HOLD | `min_frame_fraction` 0.8 |
| A5 | p2p / pipeline error below 0.05 (degenerate): REJECT; p2p above 0.05: HOLD | `min_noise_ratio` 0.05, `max_p2p` 0.05 |
| R1 | finite points below 98% of v2's: HOLD | `min_points_ratio` **0.98** (was 1.0 in PR #16) |
| R2 | p2p more than 1.2 times v2's: HOLD | `max_scatter_ratio` 1.2 |
| R3 | v2 is the plan's star and the candidate is not, or v2 has the target's aperture 5 and the candidate none: REJECT | |

The promotion step adds three:

| Check | Holds when |
|---|---|
| **S1** | a file in the candidate predates the run: its inode last changed before the run started. Reruns reuse v3 directories, so an earlier run's files can sit beside the new ones (a stale aperture-8 curve, an old `_output.fits` of another star, last run's lightcurve plots), or make up the whole directory when the latest run did not get as far as this target. It uses ctime, not mtime: T8 restores stacks from `v2/StackImages` with `shutil.copy2`, which keeps their old mtime, and on 2026-10-17 mtime would have held 16 more of the 169 v3 directories for that alone. The approval policy (below) tells clutter from stale products, and a promotion sets these files aside. |
| **S2** | no run covering the target is found, so S1 cannot be checked. |
| **A3** (with `--identity table`) | the plan file names another star than the K1 target table `ml_40pc_v2.csv` (with both of its Gaia IDs) for an Sp target. |

`--identity plan` (the default) is the checker's behaviour: plan file line 4, else the 40 pc list by name.
`--identity table` takes Sp targets' IDs from `<basedir>/ml_40pc_v2.csv`, when that file is there, and keeps the
plan file for other targets. `--target-table` points elsewhere. The table is not in the data root yet (K1, PR #17,
is not deployed), so in production `--identity table` falls back to plan files, and says so.

## PROMOTE

In this order, each step journalled before and after it runs:

1. **check**: the night is not busy (below); the candidate and the incumbent are byte for byte as they were judged
   (a fingerprint of every file's path, size and mtime, and the directory's inode); every directory involved is on
   one filesystem and writable. If this fails, nothing has moved: the promotion is *aborted*.
1a. **aside**: files in the candidate that predate its run (what S1 found) are moved, with their relative paths, to
   `superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-p<id>-stale/` with a `superseded.json`, so that v2 gets only the
   promoted run's own files. Skipped when every file is the run's own (an automatic PROMOTE always is).
2. **supersede**: the v2 incumbent, if any, is renamed to
   `PipelineOutput/superseded/<TEL>/<NIGHT>/<TARGET>/<UTC stamp>-p<promotion id>/`. That directory is created empty
   first, and creation fails if the name is taken, so a supersession can never land inside another, however many
   times a night is redone.
3. **publish**: the candidate is renamed into v2. If this fails, the incumbent is renamed back.
4. **touch**: its `*_5_diff.fits` get the current mtime, so the portal's 90-day window ingests them.
5. **reduction**:
   - When v2 has no reduction for the night, or one without masters, a copy of the v3 reduction is installed.
   - Otherwise v2's is replaced only when the run covered every target of the night (its `1_image_<target>.list`
     files), the v3 reduction has masters and image lists, and its masters were written during the run.
   - The copy is made beside v2's reduction first, then the old one is superseded under
     `superseded/<TEL>/<NIGHT>/reduction/<stamp>/` and the copy renamed into place, so v2 is never left without a
     reduction.
   - The v3 reduction stays where it is.
6. **records**: the run's log and the night's PDF (if this run wrote it) are copied, never moved or overwritten,
   to `v2/<TEL>/logs/<NIGHT>_1_v3[_<targets>].run<start>[-j<job>].log` and
   `v2/<TEL>/reports/<TEL>_<NIGHT>.run<start>[-j<job>].pdf`. When v2 has no `<TEL>_<NIGHT>.pdf` (or plain log) yet,
   that name is made a hard link to the copy, because `completeness_check.py` (05:00 cron) treats a night as
   processed only if `v2/<TEL>/reports/<TEL>_<NIGHT>.pdf` exists.
7. **thin**: the superseded copy, and the set-aside files, lose their frames and per-frame files: every file directly
   inside their run directories (`1/`, and any other all-digit one), apart from the frame lists (`*.dat`,
   `*.dat_phot`); a superseded reduction loses `flats/`. `*_diff.fits`, `*_output.fits`, `lightcurves/`, plots and stacks stay.
   `thinned.tsv` beside them lists every file deleted, and `superseded.json` says what replaced the copy. Raw is
   on disk, and v2_old was thinned the same way on 2026-10-10.

The frames go last, so v2 is without the target only for the moment between two renames rather than for the time
the deletions take.

Nothing is ever copied across filesystems or deleted outside `PipelineOutput/superseded`.

## HOLD and REJECT: the approval policy

With `--apply` the verdict, its reasons, the measurements and the fingerprints go into the journal. Then the
approval policy (`src/utils/promote_policy.py`, on unless `promote.policy` is false in the queue config or the run
has `--no-policy`) decides it, so that a person sees only what the evidence cannot settle. It answers every pending
HOLD or REJECT, new ones and ones the journal already holds, once per policy version, and journals each answer in
`policy_decisions` with its rule and evidence:

| Outcome | What happens | When |
|---|---|---|
| **promote** | the PROMOTE steps, as `accept` would; the verdict becomes *accepted*, decided by `policy:v1:<rule>` | S1 found only clutter (plots, per-frame sidecars, frame lists and frames, stacks and catalogues, an SP-named `_output.fits` of the old naming, curves of other apertures) and the checks passed; or the candidate has at least v2's points and no more scatter, or is within 5% of v2 on both; or v2 has nothing usable and the curve is a short complete night, a partial night whose missing frames the sky explains (their plate solves found too few stars), or noisy on a night with most of its frames flagged bad weather |
| **reject** | closed; v2 keeps what it has | no aperture-5 curve, no finite point, a degenerate curve, v2 the plan's star and the candidate not, v2 better on points and scatter, a few points of a long night, or scatter far beyond a light curve (p2p above 0.25) |
| **reprocess** | stays pending, listed for a rerun, not for a person (`q promote list --needs reprocess`) | the processing is at fault: frames lost to plate-solve failures on frames with plenty of stars (more than 10% of the raw frames, and at least 5), a directory the latest science run (almost) did not write, or a curve the portal serves that predates the run beside its new one |
| **person** | stays pending for a person | mixed against v2 (fewer points but less scatter, or slightly worse on both), an identity disagreement (A3), no known run (S2), a partial night the log does not explain, or a noisy curve on a clear night |

The deciding question for a partial or noisy curve is whether the data or the processing made it so: noisy is fine
when the target is right and the data are genuinely noisy, not when the pipeline processed them poorly. The run's
log says why each frame missing from the curve dropped out (astrometry.py's per-frame summary: a solve that found
too few stars is the sky; one rejected, for instance by the plate-scale check, on a frame with 20 or more sources is
the pipeline), and the curve's bad-weather flags say the same for noise. Thresholds are in
`promote_policy.THRESHOLDS` and can be overridden in `promote.policy_thresholds`.

On the 127 verdicts pending on 2026-10-19 it gave 56 promote, 20 reject, 49 reprocess and 2 person
(`mh_scratch/promote_v3_20261019/policy_dev/`). Nearly all the reprocess answers are "Plate scale sanity check
failed" on frames with hundreds of stars, from the `--force-platesolve` reruns.

A policy reject is left alone while the candidate's files are unchanged, like a person's. What waits, and the
commands for a person:

```bash
q promote list                     # what waits, and any promotion that did not finish
q promote list --needs person      # only what the policy left for a person (--needs reprocess: for a rerun)
q promote show 12                  # one verdict: reasons, notes, measurements, run, policy decisions, promotion steps
q promote accept 12 --reason "short night, fine"            # shows what it would do
q promote accept 12 --reason "short night, fine" --apply    # does it (the PROMOTE steps)
q promote reject 12 --reason "keep v2's aperture-5 curve"   # closes it; nothing moves
q promote list --grep S1                                     # only the holds for stale files
q promote accept 12-20,25 --reason "stale plots only" --apply   # several at once (one process per verdict)
```

`accept` refuses if the candidate or v2 changed after the verdict; run the night again to judge it anew. A
candidate rejected by a person, or by the policy, is left alone on later runs while its files are unchanged. A held candidate that
is judged again unchanged keeps its verdict, rather than getting a new one each run. A verdict whose candidate
changes, or leaves v3, becomes *outdated*.

## Crashes

The journal (`<queue_root>/promote.sqlite`, SQLite in WAL mode on appct's local disk, opened only inside
orchard-server) records each promotion's plan, with every path decided before the first step, and each step's
intent and outcome. Each step can be run again: the renames recognise their own result by inode, the reduction copy
is redone, and records and thinning skip what is already done.

- A promotion cut short (the process killed, the container restarted) stays *in_progress*. `run --apply` on that
  night finishes it before judging the night again; so does `q promote resume --apply`. Until then the night is
  skipped, and the daily report names the promotion.
- A step that fails after something moved leaves the promotion *failed*. The night waits for a person:
  `q promote resume --promotion N --apply` retries it, or `q promote abandon N --reason "..."` closes it once it
  has been sorted out by hand.

## Locks

A night is left alone, and a queue job for it exits 75 (retried an hour later), while any of these work on it:

- a queue job holding `night:<TEL>:<NIGHT>` (pipeline, fetch, the 16-bit frame sweep);
- a pipeline, download or promotion process started outside the queue.

The promote job holds the night lock itself. A promotion run outside the queue that changes the tree carries
`--apply --telescope T --night N` on its command line, and the dispatcher's scan of `/proc`
(`jobqueue.resources.parse_promote_argv`) counts that as holding the night. `run --apply` over several nights, and
`accept` and `resume`, therefore run one child process per night. Dry runs only read and hold nothing.

## Raw-frame counts

A4 counts the raw science frames of each night from their headers (`promote_check.RawFrames`). The counts are kept in
`<queue_root>/promote_raw.sqlite`, which concurrent promote jobs share. An entry is dropped when the night's raw
directory, or a subdirectory, has changed since it was counted. `--raw-seed` starts from `promote_check`'s JSON cache,
whose entries count as made when that file was written. The ESO backfill added about 390,000 frames to 845 nights
after the 2026-10-09 cache, so those nights are counted again. A dry run only reads the default cache.

## Running it

Inside orchard-server (`q promote ...` is the same as `python utils/promote.py ...`):

```bash
python utils/promote.py run                                   # dry run over all of v3
python utils/promote.py run --telescope Europa --night 20251016 --report dry.md --csv dry.csv
python utils/promote.py run --apply --telescope Europa --night 20251016
python utils/promote.py run --apply                           # every night in v3, one child per night
```

## In the queue

- `q add --class P2 --telescope Europa --night 20251016 --promote` queues the pipeline job with `--no_T12` and a
  promote job after it. `--promote` works with `--download` and `--from-file` too, so it is how the ESO backfill
  reruns go: `q add --class P2 --from-file backlog.csv --promote --campaign "ESO backfill reruns"`.
- The promote job (`python utils/promote.py run --apply --telescope T --night N --job <pipeline job>`):
  - runs whether the pipeline job succeeded or failed, and judges the targets that job ran (all of the night's
    if it was given none). It is cancelled if the pipeline job is cancelled.
  - holds the night lock; it is 2 cores, not disk-heavy, about 15 minutes.
  - runs ahead of the pipeline jobs of its class, because it holds up their results.
  - exits 0 when it judged everything, whatever the verdicts; 1 if a promotion failed; 75 if the night was busy.
- `q requeue` of a failed pipeline job requeues its finished promote job too, so the rerun's output is judged.
- `q add --class P2 --kind promote --after 1544-1670` adds promote jobs after pipeline jobs already queued or
  finished (with `--no_T12` or not); `q add --class P2 --kind promote --telescope T --night N` judges a night on its
  own, against its latest run.
- `promote.nightly` in `queue/config.json` (off) makes the watcher's, the look-back's and fetch reruns' pipeline jobs
  run with `--no_T12` and a promote job after them: the switch for moving the nightly runs off T12.

## Not done yet

Phase 1 promotes from v3 as it is. Per-run staging directories and provenance keywords in the products come later,
after the backfill. Also not done:

- **v3 cleanup.** T12 deleted a v3 night directory left with only its reduction. Promotion leaves v3's reduction
  and the emptied night directories alone.
- **The digest.** The 07:30 report lists holds and rejects; the morning digest (notifier) does not yet.
- **The nightly runs** still go through T12 until `promote.nightly` is turned on.
