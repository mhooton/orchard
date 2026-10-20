# Calibration register

T1 (`src/main/createlists.py`) chooses each night's bias, dark and flat frames from a register of every
calibration frame in the archive, instead of walking back one night at a time.

## Why

The old search looked for missing calibrations backwards, night by night, re-reading every header of every night
it passed: up to 3 years for flats, 100 nights for bias and darks. It took whatever it reached first, from any
camera, at any exposure time, and for flats added nights together until it had four.

- Callisto 20250407 (SPIRIT, queue job 1043) spent 2h20m in T1 and built its Clear master flat from ANDOR flats of
  20230908 (2048 x 2042) for 1024 x 1280 SPIRIT frames; the J filter, which has no flats at all, cost 1h25m.
- v2 SPIRIT nights 20220111 and 20220607 took ANDOR bias and darks (from 20220107 and 20220301). pipered.py
  sets dark = 0 when the shapes differ, so those frames were left uncorrected with only a warning.

## The register

`src/calibration/calregister.py`, one SQLite file, by default `<data root>/calibration_register/register.sqlite`
(`$CALIB_REGISTER` overrides it for T1; the queue uses its `calib_register` setting). One row per calibration
frame, keyed by telescope, night directory and the frame's path without FITS suffixes (a frame that later gains
`.fz` keeps its row):

| group | columns |
| --- | --- |
| what | `kind` (bias/dark/flat, T1's IMAGETYP rule), `filter`, `exptime`, `date_obs` |
| configuration | `naxis1`, `naxis2`, `xbin`, `ybin`, `xorg`, `yorg` (window), `orient` (MaxIm `FLIPSTAT` or Astra `ASTRAROT`/`ASTRAMIR`), `config` (all of these), `pixsize` (XPIXSZ / XBINNING), `detector` |
| detector state | `ccd_temp`, `set_temp`, `readmode`, `os_level`, `os_noise` (the overscan signature, `readout.py`) |
| level | `level` (median of 64 central rows minus the overscan), `p999`, `sat_frac`, `saturation` (raw ceiling) |
| usability | `zero` (the night's bias level), `signal`, `usable`, `reason` |

and one row per night (`nights`): file count, calibration counts, the night's detector, and a summary of a few
science frames (configuration, filters, exposures, temperature, the WCS's sky angle).

Why SQLite rather than CSV: ~600,000 rows, read by every T1 (several nights per telescope run at once since
PR #39) and updated by T1, the queue's download and fetch jobs and the builder, so it needs indexed queries and
atomic updates from concurrent writers. Like the queue's database it must only be opened from orchard-server on
appct: appcs sees that disk over NFS, where SQLite locking is unsafe. `export` writes a CSV copy.

### Usability (RULES_VERSION 1)

`flag_night()` judges each frame against the other frames of its night; `reflag` re-applies the rules without
reading any frame.

- unusable if 5% or more of its central band is at the raw ceiling (hot pixels and stars in twilight flats are far
  fewer): readout.py's saturation for the Andor readout modes; for SPIRIT 14,400 ADU, where its pixels fill (the
  configured 11,000 would call every dark longer than 30 s and most twilight flats saturated);
- bias and darks: unusable if further than max(30 ADU, 6 robust sigma) from the median of the night's frames of the
  same configuration, readout mode (and, for darks, exposure);
- flats: unusable below 1,000 ADU of signal or with a raw median above 85% of the ceiling.

## Choosing (`src/calibration/calselect.py`)

For each type, one night: the nearest to the science night on either side (the past wins a tie), among nights
with enough usable frames that are compatible with the science frames:

- same telescope and **configuration epoch**: the same camera run (a run ends wherever the night's detector -
  unbinned pixel size, else shape - changes, so Callisto's ANDOR periods either side of SPIRIT are different
  epochs) and the same side of every break in `src/calibration/calib_epoch_breaks.csv`;
- the same `config`: shape, binning, window and orientation;
- bias and darks in the science frames' readout mode; darks within the temperature tolerance;
- minimum usable frames (`calselect.py`): bias 5 (Andor nights take 6-14; with the overscan removing the level,
  a neighbour's 6-14 beat a night's own 2-4), 2 for SPIRIT (nights of 2025-26 take 1-3, and its zero level
  cancels against exposure-matched darks); darks 3; flats 4 per filter.

SPIRIT (`dark_matching_exptime`) darks are chosen per science exposure time, rounded to whole seconds as
pipered.py matches them. One night that covers every exposure that has darks anywhere is used for all of them if
it is within `COVERING_SLACK_DAYS` (3) of each exposure's own best night. The dark list also keeps that
(primary) night's darks at exposures no science frame needs, as the old lists kept all of a night's darks, so
`MasterDark.fits` - which pipeflat.py scales to every flat and pipered.py to an unmatched exposure - stays the
night's longest master; an exposure already taken from another night is never topped up, so no master mixes
nights.

Nothing is combined across nights. When nothing matches:

| case | result |
| --- | --- |
| no bias, no dark at all, no flat for any science filter | T1 exits: `ERROR: T1 cannot calibrate <tel> <night>: ...` (the night stops; the queue records the signature) |
| a science filter without flats | flagged; T6 leaves that filter's frames unreduced (as it always did) |
| a SPIRIT exposure without darks | flagged; T6 subtracts the highest-exposure master scaled to it (`MasterDark.fits`) |
| ANDOR darks only outside the temperature tolerance | used, flagged |
| science frames in a second configuration or readout mode | flagged |

Every choice is printed in the T1 log with its distance in days and written to
`reduction/<run>_calib_selection.json`; the PDF report prints one line under the master frames (red when
something was flagged).

## Keeping it current

- T1 registers its own night (pixel statistics of its calibration frames) before choosing, unless
  `CALIB_REGISTER_READONLY=1`;
- the queue's download and fetch jobs register a night that gained frames (only when the register exists; a
  failure never fails the job; `calib_register_hook: false` turns it off);
- `python -m calibration.calregister update --days 30` registers nights that are missing or whose file count
  changed (Artemis nights arrive by transfer, not through the queue's download).

## Commands

All read frames only and write only the register; run inside orchard-server.

```
python -m calibration.calregister eso-index --out <dir>/eso_index.csv     # ESO's calibration rows (prefilter)
python -m calibration.calregister build --eso-index <dir>/eso_index.csv --workers 8 --rate 100
python -m calibration.calregister update --days 30
python -m calibration.calregister night <Tel> <night>
python -m calibration.calregister validate <Tel> <night,...> --eso-index <csv>
python -m calibration.calregister reflag | summary | export <csv> | rotations <Tel>
python -m calibration.calcompare --register <reg> compare --out <dir> [--science spirit_science.csv]
```

The builder does not open every frame (11.1 M files): it reads ESO-named frames that ESO lists as BIAS, DARK or
FLAT (on SPIRIT datacube nights, frames unpacked from calibration cubes), every non-ESO frame in a subdirectory or
with a calibration name, and three science frames a night. `validate` compares this with full header scans.

## Deploying

1. Build the register at its production path (or copy a scratch build and run `update` for the nights since).
2. Deploy the code (`q release deploy`). From then on T1 refuses to run without a register.
