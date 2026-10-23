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

## Parked elsewhere

- Issue #31: monitor readout mode, gain and read noise per night and make them visible.
- The static aperture offset of ~(0.35, 0.7) px between the stack catalogue's X/Y and its WCS, on every night.
- PWV reprocessing of the ~3,541 nights whose PWV correction was missing before 2024. It was scoped at ~1,100
  CPU-hours; the code is on main, and the by-night/by-target decision is open.
