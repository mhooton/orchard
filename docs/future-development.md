# Possible future development

Things worth doing in orchard that nobody has started. Add to the list freely; when work starts, move the item into
a PR or an issue and delete it here. Started 2026-10-23, when development paused for Matthew's leave.

## Calibration

- **SPIRIT dark model.** Fit a per-pixel offset + rate to each night's dark series and evaluate it at every flat's
  and science frame's exposure, in place of scaling the longest master. Today that scaling leaves 5-11 ADU
  (0.06-0.17 % of a twilight flat) of pixel-level error, part master-bias noise and part a real ~5 ADU offset
  pattern. Most SPIRIT science exposures (1-102 s) have no matching dark at all. The measurements are in the
  write-up on PR #85 (`mh_scratch/calib_register/spirit_dark_scaling.py`).
- **Readout mode of binned and windowed Andor frames.** The register and readout.py give them mode `DEFAULT`, because
  the full-frame overscan slice does not fit. So in Artemis's 3.8 e-/ADU periods a binned bias of the other mode
  cannot be told apart. Measuring the binned overscan rows would let bias and darks match on mode there too.
- **Whole-night faults in the register.** The usability rules judge a frame against its own night, so a night whose
  every bias is off (a BZERO-state fault like Callisto March 2026) passes. A check of each night's levels against
  its epoch's would catch it.
- **Keeping the register current without T1.** Nights that arrive outside the queue's download and fetch jobs
  (Artemis transfers) are registered only when T1 runs on them, so a night not yet processed cannot lend its
  calibrations. `python -m calibration.calregister update --days 30` in the daily look-back would fix that.
- **Nights the register cannot calibrate** (T1 stops with an error):
  - Callisto 20220108, 20220115, 20220117: frames stored "Flip/Rotate 90 CW", and no calibrations exist in that
    orientation. Either turn the frames to the stored orientation (repair_orientation.py) or calibrate with
    transformed masters.
  - Artemis 20190930 and 20250704: windowed science. Masters could be cut from full-frame calibrations of the same
    binning.
  - Artemis 20210414: 2048x2048 frames without overscan.
- **SPIRIT Clear flats.** None have been taken since 2025-02, so Clear science frames are not reduced. This is
  operations (ask for Clear flats) or a policy decision.
- **Distance limits.** Calibrations more than 30 days away are used and flagged. A hard limit per type, or a
  different number, may be wanted once the flags have been watched for a while.

## Job queue

- **`q add --kind promote` ignores `--dry-run`.** `cmd_add` (`src/jobqueue/cli.py`) hands `--kind promote` to
  `cmd_add_promote` before any `--dry-run` check, and `cmd_add_promote` never reads the flag. So
  `q add --class P2 --kind promote --after 1544-1670 --dry-run` really adds a promote job after each pipeline job it
  names, and one for the night without `--after`. Found during the 2026-10-19 sweep of v3 into v2. The fix: with
  `--dry-run`, `cmd_add_promote` prints what it would add, including where a promote job with the same dedupe key
  (`promote:<tel>:<night>:<job or *>`) is already queued or running, and adds nothing. A test beside
  `test_q_add_promote_queues_a_pipeline_run_without_T12_and_a_promote_job` in
  `src/jobqueue/tests/test_promote_queue.py` should check that a dry run leaves the jobs table as it was. Until then,
  check chains with a read-only query of queue.sqlite (`mh_scratch/promote_v3_20261019/chain_check.py`).

## Telescope faults

- **ACP's pointing re-slews in the clock offset.** `MountLog.clock_offset()` (`src/utils/telescope_faults.py`)
  gives -42 to -45 s on 38 Artemis nights where the position log's clock in fact runs 1-3 s ahead of ACP's. After
  an autofocus, FocusMax slews to a focus star and back with no `Start slew to` line in ACP's log. ACP then takes a
  pointing exposure and logs `[sync] telescope synced`, `Re-slew to target.` and `Start slew to <target>...`
  (2021-09-10: acquisition 20:48:03, FocusMax's return slew 20:51:59, re-slew 20:52:46). The re-slew moves about an
  arcminute, which the log does not show as a slew, and comes 4-5 min after the acquisition, outside the 150-s
  same-target rule. So it pairs with the return slew and votes for -42 to -45 s, and on nights with two or three
  targets those votes win. The fix: skip a `Start slew to` line that follows `Re-slew to target.` (5,384 of the
  lines in Artemis's ACP logs do), keep the 150-s rule, and add a test with a synthetic Artemis-like night. Check
  also the slews after a pointing exposure that did not solve (`Start slew to offset...`, then back to the target;
  825 lines), which may mispair the same way. Validate in the container: 2021-09-10, 2021-12-13 and 2024-05-13
  should come out near +2 s; 2021-09-12, 2023-08-15 and 2024-02-11 must stay at +2.0; Europa 2019-04-19 to 04-22
  at -17 to -35 s and 2019-05-05 at +31 s. The module's offset for every Artemis night is in
  `mh_scratch/artemis_tracking/census.jsonl`, an independent long-slew estimate in `clock_check.jsonl`, and
  `slewpairs.py` there shows the pairing. With the corrected offsets the census finds the same stopped frames on
  every night, so no report is wrong today; the fix keeps T11 right on archive reruns of MaxIm nights.

## Target identification

- **Wrong-field nights go looking for a backup.** When the telescope pointed somewhere other than the target, T8
  reports `PRIMARY_NOT_IN_FOV+COORD_NO_MATCH` and then restores a backup stack catalogue of the target's real
  field, which can never land on the frames. Callisto 2024-10-28 to 10-31 (Sp0351-0052) pointed at Dec +0.879
  instead of -0.879, 1.76 deg from the star (all four nights are in `docs/unrecoverable_nights.csv`). The rerun of
  2024-10-29 (queue job 961) restored Io 2021-11-28's v2 catalogue, and condense crashed on it. Since #38 v2
  catalogues are refused and that night stops at T8, but a current-format backup would still be restored, and T9
  would be asked to measure positions 1.76 deg off the frames. The fix: before the backup search, test whether the
  target's expected position (the coordinate fallback already propagates it with proper motion) lies inside the
  solved footprint of tonight's stack. If not, stop with a status of its own (`PRIMARY_OUTSIDE_FIELD`, say),
  search for no backup, and say in the log that the night belongs in the unrecoverable list. A test can use the
  2024-10-29 stack's WCS (CRVAL 57.771, +0.875) with the target at (57.750, -0.882).
- **Catching declination sign losses when they happen.** All four Callisto nights of a -00 deg target in 2024 lost
  the sign. The other telescopes' 2024 nights of such targets are right (Ganymede 7, Io 3). Callisto has observed
  no -00 deg target since, so whether it still happens is unknown; worth raising with SSO. A check in the daily
  look-back of each night's header Dec sign against the Sp name (or the target list) would flag a new case the
  next morning. The header scan (`reports/target_catalogue_stage1.txt`, 8 Oct) has four more Callisto nights
  whose sign disagrees with the name, none checked or listed yet: Sp1500-0039 on 2018-06-02 (+0.658, another
  -00 deg target), and three whole-sign flips in Callisto's first weeks: Sp1735+2634 2018-03-26, Sp1646+3434
  2018-03-27 and Sp1235-1318 2018-04-23.

## Parked elsewhere

- Issue #31: monitor readout mode, gain and read noise per night and make them visible.
- The static aperture offset of ~(0.35, 0.7) px between the stack catalogue's X/Y and its WCS, on every night.
- PWV reprocessing of the ~3,541 nights whose PWV correction was missing before 2024. It was scoped at ~1,100
  CPU-hours; the code is on main, and the by-night/by-target decision is open.
