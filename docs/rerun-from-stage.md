# Rerun from stage N: fixing a late stage of nights already in v2

A fix to a late stage of `ZLP_pipeline.sh` should not cost a rerun from raw. When the PWV correction fails (T10's
LHATPRO fetch, as during appct's DNS outage of 2026-10-18..20), T1-T9 are fine and only T10 and T11 need to run
again. Before this, that was possible only while a night was still in v3: the promotion step moves a night's target
directories into v2, so no stage could be rerun afterwards, and a partial rerun in v3 was held by S1, because every
file it did not rewrite predates it.

The same goes for the stale-PC reprocess (forced apertures placed with a stale PC matrix: a T9 fault, fixed on
2026-10-17), planned as ~1,200 full reruns only because the frames had already moved into v2.

This adds four things, for reruns from T8, T9, T10 or T11:

1. **The stage map** (`src/utils/stages.py`): what each stage reads and writes, so that "PWV failed" resolves to
   T10 + T11 and "stale PC" to T9-T11, and a rerun's inputs and outputs are known file by file.
2. **`--from_TN`** in `ZLP_pipeline.sh`, and `q add --from TN` (or `--from pwv`, `--from stale-pc`) in the queue.
3. **A restage step** (`src/utils/restage.py`, a `restage` queue job): puts exactly the inputs of the rerun back
   into a fresh v3 target directory, as hard links to v2's own files where the rerun only reads them and as copies
   where it edits them (the frames, from T9 back), and records them.
4. **An overlay promotion**: the promotion of a restaged rerun replaces only the files the rerun wrote (its
   manifest) in v2, sets the v2 files they replace aside, and judges S1 on the manifest only.

## The stage map

`python utils/stages.py` prints it; `--cause pwv` or `--from 10` says what a rerun runs, reads and rewrites;
`--markdown` prints the tables below (a test fails if this page and the code disagree). Scopes: `raw`
Observations/&lt;TEL&gt;/images/&lt;NIGHT&gt;/, `obs` Observations/&lt;TEL&gt;/, `reduction`
v3/&lt;TEL&gt;/output/&lt;NIGHT&gt;/reduction/, `target` v3/&lt;TEL&gt;/output/&lt;NIGHT&gt;/&lt;TARGET&gt;/, `telescope`
v3/&lt;TEL&gt;/, `data` the data root. `(Tn)` after a read is the stage that made the file (T4 stands for the
masters of T2-T4). Every stage also adds rows to `reports/Pipeline_Report` and writes the night's log.

| Stage | ZLP function | Script | Reads | Writes | Edits in place |
|---|---|---|---|---|---|
| T1 lists | `create_input_lists` | `main/createlists.py` | `raw` * | `reduction` 1_bias.list<br>`reduction` 1_dark.list<br>`reduction` 1_flat_*.list<br>`reduction` 1_image_*.list | - |
| T2 master bias | `create_master_bias` | `calibration/pipebias.py` | `reduction` 1_bias.list (T1)<br>`raw` * | `reduction` 1_MasterBias.fits<br>`reduction` overscan.dat<br>`reduction` readout.json<br>`reduction` readoutnoise.dat | - |
| T3 master dark | `create_master_dark` | `calibration/pipedark.py` | `reduction` 1_dark.list (T1)<br>`reduction` 1_MasterBias.fits (T2)<br>`raw` * | `reduction` 1_MasterDark*.fits | - |
| T4 master flat | `create_master_flat` | `calibration/pipeflat.py` | `reduction` 1_flat_*.list (T1)<br>`reduction` 1_MasterBias.fits (T2)<br>`reduction` 1_MasterDark*.fits (T3)<br>`raw` * | `reduction` 1_MasterFlat_*.fits<br>`reduction` datafile.dat<br>`reduction` expdata.dat<br>`reduction` flats/** | - |
| T5 bad pixel map | `create_bad_pixel_map` | `calibration/pipebadpixel.py` | `reduction` 1_image_{target}.list (T1)<br>`reduction` 1_Master*.fits (T4) | `reduction` 1_BadPixelMap.fits<br>`reduction` hot_pixel_map.fits<br>`reduction` cold_pixel_map.fits<br>`reduction` flatbad_pixel_map.fits<br>`reduction` darkcurrent.dat | - |
| T6 calibrate | `reduce_science_images` | `calibration/pipered.py` | `reduction` 1_image_{target}.list (T1)<br>`reduction` 1_Master*.fits (T4)<br>`reduction` 1_BadPixelMap.fits (T5)<br>`reduction` *.dat (T5)<br>`reduction` readout.json (T2)<br>`raw` * | `target` 1/proc*.fits<br>`target` 1/*_processed.dat | - |
| T7 astrometry | `check_astrometry` | `astrom/astrometry.py` | `target` 1/*_processed.dat (T6)<br>`target` 1/proc*.fits (T6)<br>`reduction` 1_image_{target}.list (T1)<br>`data` Gaia database | - | `target` 1/proc*.fits<br>`raw` * |
| T8 stack | `create_stack_image` | `photometry/catalogue_fov.py` | `target` 1/*_processed.dat (T6)<br>`target` 1/proc*.fits (T6)<br>`obs` schedule/**<br>`data` ml_40pc_v2.csv<br>`data` PipelineOutput/v2/StackImages/** | `target` 1/*_stacked.dat<br>`target` 1/*.cat<br>`target` 1_initial-catalogue.fits<br>`target` *_outstack_*.fits<br>`target` *_stack_catalogue_*<br>`target` .expected_stack_catalogue<br>`target` .active_stack_catalogue_*<br>`data` PipelineOutput/v2/StackImages/** | `target` 1/proc*.fits |
| T9 photometry | `perform_aperture_photometry` | `photometry/ZLP_app_photom.py, condense/condense_photometry.py` | `target` 1/*_processed.dat (T6)<br>`target` 1/proc*.fits (T6)<br>`target` .active_stack_catalogue_* (T8)<br>`target` .expected_stack_catalogue (T8)<br>`target` *_stack_catalogue_*.fits (T8)<br>`obs` schedule/**<br>`data` ml_40pc_v2.csv<br>`data` Gaia database | `target` 1/*.cat<br>`target` 1/*.phot<br>`target` 1/*.ell<br>`target` 1/*.png<br>`target` 1/*_processed.dat_phot<br>`target` 1_initial-catalogue.fits<br>`target` *_stack_catalogue_pm_*.fits<br>`target` *_{night}_output.fits<br>`telescope` output/StackImages/already_run.txt | `target` 1/proc*.fits<br>`target` *_stack_catalogue_*.fits<br>`telescope` output/*_output.fits |
| T10 light curves | `perform_differential_photometry` | `diff_photometry/SPlightcurve.py` | `target` *_{night}_output.fits (T9)<br>`reduction` 1_image_{target}.list (T1)<br>`raw` *<br>`data` ml_40pc_v2.csv<br>`data` pwv/PWV_grids/*<br>`data` tests/targetListTeffEsts.csv<br>`data` ESO LHATPRO | `target` *_diff.fits<br>`target` *_pwv_original.fits<br>`target` *_pwv_processed.fits<br>`target` *_pwv_interpolated.fits<br>`target` lightcurves/** | - |
| T11 PDF report | `pdf_report` | `reporting/pdf_report_catriona.py` | `reduction` 1_Master*.fits (T4)<br>`target` *_output.fits (T9)<br>`target` *_diff.fits (T10)<br>`target` lightcurves/** (T10)<br>`target` 1/proc*.fits (T6)<br>`obs` technical_logs/{day}/*<br>`obs` technical_logs/log_dome_closure.txt<br>`obs` schedule/Archive_night_blocks/* | `telescope` reports/{tel}_{night}.pdf<br>`telescope` output/StackImages_nightly/{night}/*<br>`obs` technical_logs/{day}/*_PSF_*.png | - |
| T12 sweep | `perform_version_migration` | `main/ZLP_pipeline.sh` | - | - | - |

| Cause | Rerun from | Stages | What |
|---|---|---|---|
| `pwv` | T10 | T10 T11 | the PWV correction: the LHATPRO fetch, the PWV grids, the correction itself |
| `teff` | T10 | T10 T11 | the target Teff the PWV correction needs (header TEFF, the Filippazzo list) |
| `comparison-stars` | T10 | T10 T11 | comparison-star weights and the artificial light curve |
| `light-curves` | T10 | T10 T11 | anything else in SPlightcurve.py: apertures 3-8, MCMC files, plots, best aperture |
| `pdf` | T11 | T11 | the night report only |
| `stale-pc` | T9 | T9 T10 T11 | forced apertures placed with a stale PC matrix over the fitted CD (fixed 2026-10-17, PR #30) |
| `photometry` | T9 | T9 T10 T11 | aperture photometry, condensing, the night's *_output.fits |
| `stack` | T8 | T8 T9 T10 T11 | the stack, its catalogue, the Gaia crossmatch, which stars are targets |
| `astrometry` | T7 | T7 T8 T9 T10 T11 | plate solving of the calibrated frames |
| `calibration` | T6 | T6 T7 T8 T9 T10 T11 | calibration of the science frames (pipered.py) |
| `masters` | T2 | T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 | master bias, dark or flat, the bad pixel map |
| `lists` | T1 | T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 | which frame is what (createlists.py) |

| Rerun from T8: input | Made by | Restaged as |
|---|---|---|
| `target` 1/*_processed.dat | T6 | copy, each path pointed at the restaged frame |
| `target` 1/proc*.fits | T6 | copy |
| `reduction` 1_image_{target}.list | T1 | hard link |
| `reduction` 1_Master*.fits | T4 | hard link |

| Rerun from T9: input | Made by | Restaged as |
|---|---|---|
| `target` 1/*_processed.dat | T6 | copy, each path pointed at the restaged frame |
| `target` 1/proc*.fits | T6 | copy |
| `target` .active_stack_catalogue_* | T8 | hard link |
| `target` .expected_stack_catalogue | T8 | hard link |
| `target` *_stack_catalogue_*.fits | T8 | copy |
| `reduction` 1_image_{target}.list | T1 | hard link |
| `reduction` 1_Master*.fits | T4 | hard link |

| Rerun from T10: input | Made by | Restaged as |
|---|---|---|
| `target` *_{night}_output.fits | T9 | hard link |
| `reduction` 1_image_{target}.list | T1 | hard link |
| `reduction` 1_Master*.fits | T4 | hard link |
| `target` *_output.fits | T9 | hard link |
| `target` 1/proc*.fits | T6 | hard link |

| Rerun from T11: input | Made by | Restaged as |
|---|---|---|
| `reduction` 1_Master*.fits | T4 | hard link |
| `target` *_output.fits | T9 | hard link |
| `target` *_diff.fits | T10 | hard link |
| `target` lightcurves/** | T10 | hard link |
| `target` 1/proc*.fits | T6 | hard link |

**From T9 and T8 the frames are copies.** Those reruns edit the calibrated frames in place: `wcsfit` (T8's stacked
frames, T9's every frame) rewrites their headers, T9 removes a stale PC matrix and sets `wcscompl`. A hard link
shares one file between v2 and the rerun, so an in-place edit would change v2's frame too; the rerun gets its own
copy of every frame (1-19 GB per target-night), and its overlay replaces the frames and their per-frame files. The
stack catalogue is a copy as well: condense opens it read-write.

**T9 needs a catalogue with TARGET_ROLE.** Condense reads the targets T8 identified from the catalogue's
TARGET_ROLE column, which T8 writes since the stack/catalogue overhaul of 2026-03-31. A night whose catalogue is older
cannot be rerun from T9 with today's code: restage says so, and the night is rerun from T8, which makes a new stack
and catalogue. Such nights also have target-named catalogues (`SP1054-5431_stack_catalogue_I+z.fits`) and no record
of the one T9 used; for a rerun from T9 of a night that has exactly one catalogue per filter, restage writes the
`.active_stack_catalogue_<filter>` record T8 writes today, naming it, and refuses to guess between several.

**Why not from T7 back.** In the code this was built on, T7 writes each frame's new WCS into the raw frame its image
list names, by absolute path, so a rerun anywhere would write production raw frames. PR #64 (on main since 2026-10-20)
replaces that with a WCS cache, which makes T7 the next stage to add; it has not been tried.
From T6 back the frames themselves are made anew, so a rerun costs what a run from raw does.

## `--from_TN` and the log

`ZLP_pipeline.sh --from_T10` runs T10 and every stage after it, less any `--no_T`; it cannot be combined with
`--only_T`. Every run now logs the stages it runs (`Stages run: T10 T11`), so the promotion step can tell, from a
log as well as from a queue job's argv, which files the run had to write. `utils/stages.stages_from_argv()` gives
the same answer as the script; a test runs the script's own argument handling against it.

Two guards for runs that skip earlier stages:

- the target-name helpers (`correct_target_names.py`, `already_run_targs.py`, which rewrites the telescope's
  `already_run.txt`) run only with T9, the only stage that uses their result;
- frame packing (16-bit) is skipped after a run from T10 or later: T10 and T11 only read the frames, and a restaged
  frame is v2's own file. A run from T9 or T8 packs its own copies, as any run does: float32 frames in v2 come back
  16-bit.

## Restaging

`python utils/restage.py stage --telescope T --night N --from 10 [--target X] --apply`, or the `restage` job that
`q add --from` queues. For each target of the night (default: every target with the rerun's inputs in v3 or v2):

- **in v2** (promoted): a fresh `v3/<TEL>/output/<NIGHT>/<TARGET>/` gets the inputs from the table above, as hard
  links to v2's files where the rerun only reads them and copies where it edits them (v2, v3 and `superseded/` are
  all on sda1: device 2049 on appct, checked 2026-10-20; across filesystems a link would be a copy). A file a stage of
  the rerun writes is never an input, whatever read pattern it matches (T9's proper-motion catalogue matches T9's
  read of T8's catalogues). The frame lists (`1/*_processed.dat`, absolute paths) are copied with each line pointed
  at the restaged frame: in production the same paths, in a scratch tree new ones. The night's v3 reduction keeps
  what it has (promotion leaves it in v3; for a single-target run it is the reduction that run made) and gets the
  image list and masters from v2 where it lacks them. v2 is only read.
- **in v3** (never promoted, or held): the directory's earlier outputs of the rerun's stages are moved to
  `superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-prerun/`, so nothing of the earlier run can mix with the rerun's (a
  stale `*_bestap.txt` that the new T10 would otherwise read and compare against, an aperture the new run fails to
  make). What stays are the earlier stages' products, the inputs.
- a v3 directory without the inputs, or a target with none anywhere, is refused.

The **record**, `v3/<TEL>/output/<NIGHT>/.restage/<TARGET>.json`, lists every input: where it came from, how it was
placed (link, copy, in place, present), its inode, size and mtime, and the v2 file's inode. It is written first
(`staging`) and marked `staged` at the end, so a restage cut short is finished by running it again, and one already
done is checked and left alone. `restage.py show [--manifest]` shows a night's records and what the rerun has
written so far; `restage.py clear --reason ... --apply` retires a restage that will not be promoted (a v2 one's
directory goes to `superseded/.../<stamp>-restage/` with its links to v2 dropped and its frame copies thinned).
A restage with `--apply
--telescope T --night N` on its command line holds the night against the queue, like a promotion.

**Why hard links are safe here.** A linked input is one the stage map says the rerun only reads, and the writers
of those stages replace a file rather than write into it (astropy's `writeto(overwrite=True)` removes the old file
first, framepack writes a temporary file and renames it). The promotion checks it anyway (O1 below): a linked input
that changed during the rerun means v2's file changed with it, and nothing moves until a person has looked. A copied
input the rerun edits by design (a frame) is part of what it wrote.

## The manifest

What a restaged rerun wrote: every file in the restaged directory that is not an unchanged input
(`restage.manifest()`), so the frames a rerun from T9 rewrote too. Each entry has its size, mtime, ctime and sha1, and for FITS files the provenance orchard
stamps on what it makes (ORCHVER and ORCHREL at creation, ORCHMOD on an edit: PR #56 and #60). The promotion journal
keeps the manifest in the promotion's plan, and `<stamp>-overlay/superseded.json` keeps it beside the files it
replaced.

## The overlay promotion

`utils/promote.py` (the promote job after the rerun) sees the record. A candidate restaged from v2 is promoted as an
**overlay**:

| Step | What |
|---|---|
| check | the night is not busy; the candidate and v2's directory are as judged; every restaged input and the v2 file it came from are as restaged; every manifest file is as listed; one filesystem; writable |
| install | each v2 file the rerun replaces, and each of v2's outputs of the rerun's stages that the rerun did not write again (an aperture it no longer makes), gets a hard link under `superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-overlay/`; then each manifest file is renamed over v2's (atomic on one filesystem, so v2 never lacks the file) and read back against its sha1; the outputs no longer written leave v2 (their data stays under `-overlay/`) |
| touch | the `*_5_diff.fits` get the current mtime, for the portal's 90-day window |
| records | the run's log and PDF are copied into v2 as `.run<start>[-j<job>]`; the PDF also takes the plain name `<TEL>_<NIGHT>.pdf` when the rerun was given every target of the night that has products (else the night's PDF would lose the others); the old plain file is kept under its own `.run` name or set aside |
| leftover | the restaged directory, now only inputs, goes to `superseded/.../<stamp>-restage/` with its record; its hard links to v2's files are dropped there |
| thin | as for any promotion; after a rerun from T9 or T8 also the replaced frames and per-frame files under `-overlay/1/` (raw is on disk, as for any superseded copy) |

v2's directory keeps its inode, and every v2 file outside the manifest is untouched: same inode, size and mtime.
The reduction is never changed (the rerun only read it). Nothing is deleted outside `superseded/`.

A night restaged in place in v3 is promoted whole, as before, with the record moved to the superseded tree after.

**S1 on the manifest.** The files the run had to write are judged, not its inputs: a restaged candidate's unchanged
inputs are exempt, and so, after a run that started at a later stage without a restage (`--from_T10`, or the
`--no_T1 ... --no_T9` reruns queued on 2026-10-20 as jobs 2064-2077), are the products of the stages it did not run.
An output of the rerun's own stages that predates it still holds, as before.

**New holds**, decided by the approval policy:

| Check | Holds when | Policy |
|---|---|---|
| O1 | an input the rerun should only have read changed or went after the restage (for a hard link: v2's file too), the record never reached `staged`, or the night's v3 reduction changed | person (`restage-inputs`) |
| O2 | a v2 file the rerun read is no longer the one restaged (v2 changed since, e.g. another promotion) | reprocess (`v2-changed-since-restage`): clear, restage and rerun |
| O3 | a FITS file the rerun made carries another code release (ORCHREL) than the run's queue job (a restaged frame it edited keeps its maker's) | person (`restage-inputs`) |

A restaged candidate is not judged until a run covering it has started after the restage (a restage from T11 holds
T10's curves from the start); one restaged but never rerun is named in the promote job's output.

## In the queue

```bash
q add --class P3 --telescope Callisto --night 20250317 --from pwv --promote --campaign "PWV outage reruns"
q add --class P3 --telescope Callisto --night 20250317 --from T10 --promote --dry-run     # says what it would do
q add --class P3 --from-file pwv_v2_nights.txt --from T10 --promote --campaign "PWV outage reruns"
q add --class P3 --telescope Io --night 20260223 --from stale-pc --promote                # T9; T8 if restage says so
```

`--from` takes a stage (T8-T11) or a cause from the table above. It queues three jobs:

1. **restage** (`utils/restage.py stage --apply ...`): 1 core, night lock; from T10 on links only (seconds, not
   disk-heavy), from T9 and T8 copies of the frames (disk-heavy, a few seconds per GB);
2. **pipeline**: `ZLP_pipeline.sh --from_TN --no_T12 ... <targets>`, given the restaged targets explicitly, so
   that its log and PDF say which targets it covers; it runs only if the restage succeeded, and has its own dedupe
   key, so it never merges with a full run of the night;
3. **promote**, as for any `--promote` run.

`--from` needs `--promote` (or `--no-T12` to keep the rerun in v3): T12 would sweep the restaged inputs into v2.
It refuses `--download`, and a `--from-file` line with DOWNLOAD=1. A rerun is queued with a share of a full run's
estimate (`jobspec.RERUN_SHARE`); at launch the dispatcher re-estimates it from the frames like any pipeline job, so
its timeout stays a full run's. Reruns from T10 or T11 reserve 8 cores (SPlightcurve.py is one process); from T9 or
T8 a full run's, since their WCS fits and photometry use every core.

## Verified

Scratch trial, 2026-10-20 (`mh_scratch/rerun_from_stage/` on appct; nothing written outside it). The scratch root
holds copies of five PWV-lost v2 nights (and what production v3 holds for them), raw frames and schedules by symlink;
the code is this branch as a `git archive` snapshot. Per night: `restage.py stage --apply`, `ZLP_pipeline.sh
--from_T10 --no_T12`, `promote.py run --apply`, and v2 compared file by file before and after (`trial.sh`,
`compare.py`; each night's report is `trial/<TEL>_<NIGHT>/compare.txt`).

| Night | Target | Frames | Stars | Restage | Rerun (T10 + T11) | Promote | v2 files untouched | Installed | Night's PDF |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| Europa 20260703 | Sp2338-4726 | 99 | 150 | 1 s | 144 s (105 + 39) | 2 s | 812 of 900 | 88 + 17 new | replaced: the only target |
| Io 20241211 | Sp0409-0605 | 135 | 146 | 1 s | 157 s (105 + 52) | 3 s | 9,199 of 9,287 | 88 + 17 new | kept: Sp0642-6707 (16.8 GB, untouched) has products |
| Callisto 20241209 | Sp0314+2840 | 379 | 53 | 2 s | 218 s (132 + 86) | 3 s | 2,509 of 2,597 | 88 + 17 new | replaced: Sp0903+0540 has none |
| Callisto 20250317 | Sp0830+0947 | 1,691 | 22 | 3 s | 478 s (367 + 111) | 6 s | 13,515 of 13,603 | 88 + 17 new | replaced |
| Europa 20250820 | Sp1804-1616 | 54 | 6,606 | 1 s | 673 s (630 + 43) | 2 s | 226 of 314 | 88 + 17 new | replaced |

- **Only the manifest changed.** On every night the files that changed in v2 are exactly the ones the promotion
  journal lists as installed (the 6 curves, 6 MCMC files, the plots and best-aperture files, and 17 new: the three
  `*_pwv_*.fits`, the `*_pwv.png` plots, `*_bestap_pwv.txt`). Every other file has its inode, size, mtime and (frames
  aside) sha1 as before; v2's directory keeps its inode; the frames' link counts are back to 1; nothing left v2. The
  replaced files are under `superseded/.../<stamp>-overlay/`, the restaged inputs under `<stamp>-restage/`.
- **The curves gained PWV and nothing else of substance.** In every replaced `_diff.fits`, IMAGELIST, CATALOGUE,
  LIGHTCURVE_n, ERROR_n and ALC_n are identical to v2's; the PWV table went from 0 rows to 453-589 and
  PWV_LIGHTCURVE_n from NaN to finite. In every MCMC file every column is identical except DIFF_FLUX_PWV, finite at
  every point now (0 before). Two more differences come from code released after the original runs (5c79fa0,
  b383cd7 and f854bc4, all before 2026-10-19 12:26 UTC), not from the rerun: a SMEARED / SMEAR_FLAG column (#55,
  the smear-trail flag) and ORCHVER in the header (#56, versioning).
- **Verdicts**: three PROMOTE; two HOLD that the policy published, as v2 itself was: Callisto 20241209 (A4, 287 of
  379 raw frames, the same points and scatter as v2: `better-than-v2`) and Europa 20250820 (15 points on a
  bad-weather night: `few-points`).
- **Dry run over all 95 production nights** (`dryrun95/`, nothing written): all 123 target-nights restage from v2,
  no refusals, 124,923 hard links; every night's v3 reduction already has the image lists and masters. The
  `--from-file` list is `dryrun95/pwv_v2_nights.txt`.
- The suite in orchard-server (67df8e9, `--basetemp` on /data): 770 passed, 2 skipped.

## What it costs: the 95 PWV-lost nights in v2

Their original full runs (the queue jobs that lost PWV, `cost/v2_nights.csv`) took **68.9 job-hours** (median 33
min a night, 20 cores each: 1,378 reserved core-hours), read **0.70 TB** of raw frames (157,229 files) and wrote
**0.46 TB** of calibrated frames, which a full rerun would do again.

A rerun from T10, from the five trial nights: each took 10-32% of its night's original run. T10 scales with the
curve's points and the field's stars, T11 with the frames; fitted to the trials (`estimate.py`), T10 per aperture =
12.3 s + 0.029 s per point + 0.014 s per star, and T11 = 47 s + 0.040 s per frame. Over the 95 nights' 125 curves:

| | Full reruns (as originally run) | Reruns from T10 |
|---|---:|---:|
| job time | 68.9 h | **~15 h** (T10 12.6 h, T11 2.6 h): 22% |
| cores reserved | 20 per job: 1,378 core-hours | 8 per job: ~120 core-hours (9%) |
| longest night | Callisto 20251224: 4 h 14 min | ~1 h |
| raw frames read | 0.70 TB | one header per target |
| calibrated frames written | 0.46 TB | none (hard links); ~105 small files per curve |
| restage and promotion | - | seconds a night (1-3 s and 2-6 s) |
| ESO LHATPRO queries | ~6 per curve | the same |

With the queue running several nights at once (they hold only their night locks), that is a few hours of wall clock.
The largest uncertainty is the crowded bulge fields of Sp1804-1616 (up to 11,072 stars; the trial's had 6,606):
the model puts their T10 at up to ~17 minutes each.

## Not done

- **Reruns from T7 or earlier** of promoted nights (see "Why not from T7 back"); T7 can follow now that PR #64 is on
  main, after a trial.
- **The archive PWV campaign** (the parked reprocessing of 3,541 nights): its nights were reduced by older
  pipelines, and their image lists can name paths that no longer exist (`/import/appct/...`); T10 opens the list's
  first frame. Restage would have to copy such a list with its paths rewritten. Their v2 directories may also lack
  frames (T11's stack panel then says so; the curves are unaffected).
- **The portal's own tables**: the overlay touches the `*_5_diff.fits` as every promotion does; whether the portal
  re-ingests a target-night it already holds is the portal's concern, as for any promotion that replaces one.
