# Telescope faults: frames the telescope lost, not the pipeline

Some nights' data are poor because the telescope lost them: the mount stopped tracking and the stars trailed, or the
pointing walked off the target. Reprocessing cannot recover those frames, and users should be able to see that it was
the telescope, not the pipeline. `src/utils/telescope_faults.py` finds both faults per target-night from data the
pipeline already has, and the night report (T11) and the morning digest show them beside the weather flags.

It reads, and never writes, the raw frames' headers, the calibrated frames' WCS (T7's solutions), T7's WCS cache
(`docs/wcs-cache.md`) and T7's per-frame lines in the night's logs. It uses the standard library only, so the job
queue can import it.

## Tracking stopped

Astra writes three keywords into every frame: the share of the exposure the mount spent tracking (`TRACKING`), the
share it spent slewing (`SLEWING`), and whether it was parked (`TELPARK`). A light frame for which the mount did
neither for at least half the exposure (`1 - TRACKING - SLEWING >= 0.5`, or parked) was exposed on a stopped mount.
Its stars trail by 15 arcsec per second.

The rest of the header agrees. `ALTITUDE`, `AZIMUTH` and `AIRMASS` stay fixed to 1e-6 or better. Before Astra 1.0 the
header RA is the mount's position, so it runs at the sidereal rate (15.04 deg/h) and Astra's own WCS is empty
(`CRVAL` 0). From Astra 1.0 (2026) RA is the target's and does not move, and Astra's WCS is frozen. On the two 2025
nights T7 failed 569 of the 570 and 332 of the 333 stopped frames.

| Measured | Stopped | Frames | Header |
|---|---|---|---|
| Ganymede 2025-07-16 TOI-4336 c | 00:30–03:11 UT | 570 of 916 | RA 15.02 deg/h, airmass 1.1187153 |
| Ganymede 2025-08-25 Sp1954+1801 | 00:32–04:42 UT | 333 of 446 | RA 15.03 deg/h, airmass 1.5547; next target 04:52 |
| Ganymede 2026-06-30 Sp2253-4801 (Astra 1.0.1) | 07:12–09:43 UT | 103 | RA fixed; altitude 58.5196, azimuth 144.5753, airmass 1.17191 |

The first frame of a block often starts while the mount still slews (`TRACKING` 0.74 with `SLEWING` 0.26). That is
not a stop. Sky flats move between exposures and are not light frames.

Astra started on Callisto on 2023-12-15, on Ganymede on 2025-02-27, on Io on 2026-01-28 and on Europa on 2026-01-29
(the first light frames with `TRACKING`; none returned to MaxIm DL). MaxIm DL frames (Artemis, and SSO before
Astra) have no `TRACKING`. Their `ALTITUDE` is computed from the target's coordinates, so their headers cannot show
a stop. For them only T7 can tell. Since v4.1.0 (`astrom/pointer_wcs.py`) T7 fails a frame whose brightest sources
are stretched one way as `Sources trailed (...)`. Three such frames in a row make a *trailed* episode.

## Pointing drift

Each solved frame's WCS puts the block's pointing (its header target position) on the detector. Guided, it stays
within a few pixels: in the archive, 98% of blocks keep 95% of their frames within 0.25 arcmin. When the pointing
walks away, the target slides across the field and can leave it.

Measured on Artemis 2024-02-11 Sp1356+4342, the pointing sat at one place to 01:33 UT, jumped 4.3 arcmin at 01:34 and
stayed there. No frame solved after 01:49 UT. The known-scale solver's re-solve of those frames had it 6.5 arcmin
off at 01:58 and about 15 arcmin off by 04:11, outside the 12-arcmin field. On Io in 2017–18 the telescope
sometimes pointed 4–17 arcmin away for a few minutes, with the target off the detector, and then came back.

A drift is the pointing more than **1 arcmin** from *home* on **3 solved frames in a row**, judged on a running median
of three so that one wrong solution is not a drift. Home is where the pointing sat longest: the median of the
densest 0.5-arcmin cell. The first frames are not used for home, because the telescope often settles first. On Io
2017-12-13 the first frames were 4.3 arcmin off and the telescope then stayed centred all night; measured from those
frames, the whole night looked displaced. Unsolved frames among the drifting ones count with them. A drift with no
solved frame after it runs to the end of the block and says so. A drift that is over within 10 minutes of the block's
first frame is the telescope settling, not a fault.

Solutions are filtered before they count:

- **Orientation.** A solution turned more than 3 degrees from the block's orientation is wrong. Old twirl solutions
  accepted before the plate-scale check land at random angles; Artemis 2024-02-11 carries ten of them.
- **Scale.** A solution whose plate scale is more than 3% off the block's is wrong.
- **Meridian flip.** The field turns 180 degrees and the pointing lands elsewhere on the detector, all at once. Each
  side of a flip has its own home.
- **Coherence.** A pointing that moves does so coherently. Guided, consecutive frames agree to a few pixels, and a
  jump is one step followed by steady frames. Wrong solutions scatter: in Callisto's 2024 products, frames 17 s apart
  land hundreds of pixels apart, with plate scales from 0.336 to 0.354 arcsec/px, inside the orientation and scale
  filters. A block whose median frame-to-frame step exceeds 0.25 arcmin cannot be judged: its drift is `unknown`, and
  its frames are counted as `noisy`. A drift episode must be as coherent itself.

An episode is a *jump* when the displacement arrives between solved frames under 5 minutes apart, and a *drift*
otherwise. The report also gives the largest displacement, when it was reached, and when the target left the field.

Each frame's WCS comes from the first of these that has one:

1. its calibrated frame: T7's solution for the run
2. T7's WCS cache row (`PipelineOutput/wcs_cache/<TEL>/<night>.ecsv`, read without astropy)
3. the raw header: Astra's own solution, or orchard's from before v4.1.1
4. T7's log line (`CRPIX`, `CRVAL`), with the matrix of the nearest frame that has a full WCS

## Where users see it

- **The night report (T11, page 2).** The *Special Events* table has two more cells beside the weather alert:
  *Tracking lost* and *Pointing drift*, `Y` in bold red, `N` in green, as the weather flags are, and `unknown`
  where the frames cannot say (MaxIm DL headers for tracking; no solved frame for drift). Each episode gets a line
  under the table, for example:

  > TOI-4336 c: tracking stopped 00:30–03:11 UT, frames 347–916 of 916 (570 frames); header RA ran 15.0 deg/h,
  > airmass fixed at 1.119

  > Sp1356+4342: pointing jumped 4.3' at 01:34 UT; 01:34–06:54 UT, frames 90–648 of 648 (559 frames); no frame
  > solved after 01:49 UT

  The same lines go into the pipeline log as `Telescope fault: ...`.
- **The record.** T11 leaves the night's result in `PipelineOutput/v3/<TEL>/reports/telescope_faults/<night>.json`.
  The name must not contain `<TEL>_<night>`: T12 moves the first file under `reports/` that matches that pattern,
  and that must be the PDF.
- **The morning digest.** The telescope's line for the night ends with what the record says, for example
  `· telescope: tracking stopped on TOI-4336 c 00:30–03:11 UT (570 frames)`. It does not count as needing
  attention: the data are poor, but nothing in the pipeline needs fixing. The digest reads the record and never
  works it out itself, because a SPIRIT night's raw headers run to 20,000 files.

## The census

CENSUS

## Thresholds

| | | |
|---|---|---|
| `IDLE_MIN` | 0.5 | share of an exposure with the mount neither tracking nor slewing |
| `TRAIL_MIN_FRAMES` | 3 | T7 `Sources trailed` frames in a row (one can be a bright star's bleed trail) |
| `DRIFT_ARCMIN` | 1.0 | pointing this far from home |
| `DRIFT_MIN_SOLVED` | 3 | solved frames in a row beyond it |
| `REF_RADIUS_ARCMIN` | 0.5 | home: the densest cell of this size |
| `SETTLE_MINUTES` | 10 | an episode over this soon after a block starts is settling |
| `ROT_TOL`, `SCALE_TOL` | 3 deg, 3% | a solution further off the block's orientation or scale is wrong |
| `JUMP_MINUTES` | 5 | a displacement arriving between solved frames this close is a jump |
| `NOISY_STEP_ARCMIN` | 0.25 | median frame-to-frame step above which the solutions, not the pointing, are moving |
| `PAUSE_MINUTES`, `PAUSE_CADENCES` | 2 min, 4 | frames this far apart mark a pause, after which settling is allowed |
| `LOW_ALTITUDE` | 25 deg | a stop below this has its altitude given (it may be the mount's limit) |
