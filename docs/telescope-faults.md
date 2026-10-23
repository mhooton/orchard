# Telescope faults: frames the telescope lost, not the pipeline

Some nights' data are poor because the telescope lost them: the mount stopped tracking and the stars trailed, or the
pointing walked off the target. Reprocessing cannot recover those frames, and users should be able to see that it was
the telescope, not the pipeline. `src/utils/telescope_faults.py` finds both faults per target-night from data the
pipeline already has, and the night report (T11) and the morning digest show them beside the weather flags.

It reads, and never writes, the raw frames' headers, the mount's own position log in the night's technical logs, the
calibrated frames' WCS (T7's solutions), T7's WCS cache (`docs/wcs-cache.md`) and T7's per-frame lines in the
night's logs. It uses the standard library only, so the job queue can import it.

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
Astra) have no `TRACKING`. Their `ALTITUDE` and `AIRMASS` are computed from the target's coordinates, so their
headers cannot show a stop. The mount's own position log can (next section). Where there is no log, only T7 can
tell. Since v4.1.0 (`astrom/pointer_wcs.py`) T7 fails a frame whose brightest sources are stretched one way as
`Sources trailed (...)`. Three such frames in a row make a *trailed* episode.

### The mount's own position log

Each telescope's ASTELCO controller writes the mount's position once a second into the night's site logs,
`Observations/<TEL>/technical_logs/<YYYY-MM-DD>/<TEL>_pos.csv`, noon to noon UT, with the columns `jd, date, time,
ALT, AZ, DEC, RA, DEC(inst.), HA(inst.), Target_Dist, Dome_AZ` (`Target_Dist` from mid-2019). The files come from
the site every afternoon with the rest of the night's technical logs. They exist for the four SSO telescopes from
February 2019 (Callisto from 02-13, the others from 02-20) and for Artemis from 2021-08-25 to 2025-06-10. Artemis's
files hold only their header from 2022-06-24 to 2023-07-28 and on 2022-03-15 to 03-18. Before 2019 the SSO site
logs hold only ESO's seeing and water-vapour logs.

Tracking, the mount's hour-angle axis turns at the sidereal rate, 0.25 degrees a minute, whatever the target's
declination, and the altitude and azimuth follow the target. When the mount stops, every axis stands still. The
altitude, the azimuth and both instrumental axes stay fixed to 1e-4 degrees or better, while the log's RA, worked out
from the fixed hour angle, runs at 15 degrees an hour. Between two lines the mount stood still if every axis moved
at less than a tenth of the sidereal rate (`STILL_RATE`). Each MaxIm DL light frame gets the share of its exposure
the mount stood still, and a frame still for half its exposure or more is a stopped frame, as with Astra's keywords.
A gap of more than 10 s between lines is unknown, and so is a frame the log covers less than half of. Where a frame
has Astra's keywords, they decide and the log is not read.

The log shows each measured stop from its first second:

| | Astra's `TRACKING` | The position log |
|---|---|---|
| Ganymede 2025-07-16 TOI-4336 c | 00:30:03–03:11:57 UT, 570 frames | still from 00:29:58 to 03:13:20 at altitude 63.32, azimuth 227.30; the same 570 frames |
| Ganymede 2025-08-25 Sp1954+1801 | 00:32:41–04:42:11 UT, 333 frames | still from 00:32:11 at 39.93, 33.98; 334 frames |
| Ganymede 2025-09-04 Sp1954+1801 | 02:59:44–04:38:58 UT, 133 frames | still from 02:59:18 at 43.07, 333.66; 134 frames |
| Callisto 2024-10-12 Sp0002+0115 | 02:04:24–04:22:39 UT, 193 frames | slewed out of the park at 02:04 and never tracked; still to 04:27:01 at 58.64, 36.55; the same 193 frames |

Frame by frame over those four nights, the log calls stopped all 1,229 frames that `TRACKING` calls stopped, and none
of the 2,083 it calls tracked. Two more are the first frame of a stop, during which the mount stood still for 60%
and 76% of the exposure. Astra's keyword still has them tracked.

**The two clocks.** The controller's clock is not the one that stamps the frames. Before Astra, ACP's PC stamped
them, and the two clocks mostly agreed to 2 s, the controller's ahead. For weeks at a time, though, they drifted
apart. Europa's position log ran 17 s behind ACP's clock on 2019-04-19 and 35 s behind by 04-22, and it ran 31 s
ahead on 05-05. Half a minute moves a stop's edge by a frame. It also puts the night's last frame onto a still
mount, though the frame ends seconds before ACP turns tracking off. On Europa 2019-04-19 to 04-22 those last frames
have round stars.

So the detector sets the log's clock against ACP's, from the slews both record: ACP's `Start slew to` lines
(`<TEL>_ACP.txt`, beside the position log) against the starts of fast motion in the position log (an axis over
0.05 deg/s). The offset most of ACP's slews agree on to 3 s is used, if two or more agree and no other offset has as
many. ACP slews to a target a second time a minute after acquiring it, to correct its pointing. The log rarely
registers so small a move, so a slew to the target of the slew before, within 150 s, has no vote. Without ACP's log
(the Astra era) or with too few slews, the offset is 0.

Measured on one night a month per telescope, from March 2019 to Astra (295 nights):

| Clock offset | Nights |
|---|---|
| log 0 to 3 s ahead | 205 |
| further apart: Europa 2019-05-15 (+36.6 s) and 2025-03-15 (−21.8 s) | 2 |
| too few distinct slews, or a tie | 25 |
| no slews in ACP's log (nothing observed) | 63 |

The nights the census read again found two more such spells: Io 2019-06-02 (−11 s), and Europa 2025-03-14 to
03-19, drifting from −16 s to −45 s at 6 s a day. Corrected, the false end-of-night stops on those nights go.

Two other checks:

- **The controller's status log.** `<TEL>_TSI_status.csv` records `CURRENT.TRACK`, the number of the trajectory the
  mount is running, or 0 when it runs none. For all 18 stops the position log found in the Astra era, and 31 of the
  79 before it, `CURRENT.TRACK` was 0 from the first frame to the last, or fell to 0 during the first frame. The other
  way it proves nothing. When an axis drive failed (Ganymede February 2022, Callisto 2021, both with `ERR_DeviceError`
  on axis 0), the controller kept a trajectory running while the axes stood still, and the frames trailed all the same.
- **The frames.** A tracked frame's stars are a few pixels across. On a stopped mount each star is a streak, 15 arcsec
  per second times the cosine of the declination. Checked frame by frame (the first, middle and last frame of each
  stop, and the frames either side; `trailcheck.py`), 62 of the 65 pre-Astra stops with enough stars detected have
  trailed frames inside, with axis ratios of 4 to 220 and up to 150 streaks; the other three show streaks by eye. The
  frame after a stop is round in 32 of 35; the frame before is round in 49 of 64, and the rest are partly trailed,
  the mount having stopped during them for less than half the exposure.

The other site logs do not show a stop on their own:

| File | During a stop |
|---|---|
| `_TSI_TrackingKernel.txt` | `Unexpected end of trajectory execution` as the three Ganymede stops began; nothing for Callisto 2024-10-12, whose mount never started tracking after its slew |
| `_events.csv` | `Data from telescope is not valid` as the three Ganymede stops began; before Astra, the causes below (`ERR_TrajectoryNoPoints`, `Time not valid. tracking points deleted`, `ERR_DeviceError` on axis 0), some of which also come on nights without a stop; `Tracking started`/`stopped` events come for some trajectories only |
| `_offsets.csv` | DONUTS stops correcting, as it also does in cloud, during an autofocus and on unguided blocks |
| `_TSI_commands.csv` | nothing until the next slew |
| `_ACP.txt` (MaxIm DL era) | the run's end (`RUN STOPPED AT SPECIFIED QUIT TIME`, `[**HARDWARE OR DRIVER ERROR**]`), not the stop |

## Pointing drift

Each solved frame's WCS puts the block's pointing (its header target position) on the detector. Guided, it stays
within a few pixels: in the archive, 97% of blocks keep 95% of their frames within 0.25 arcmin. When the pointing
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
  where the frames cannot say (MaxIm DL headers without the mount's position log for tracking; no solved frame for
  drift). Each episode gets a line under the table, for example:

  > TOI-4336 c: tracking stopped 00:30–03:11 UT, frames 347–916 of 916 (570 frames); header RA ran 15.0 deg/h,
  > airmass fixed at 1.119

  > Sp1223-2757: tracking stopped 00:44–04:27 UT, frames 122–417 of 417 (296 frames); the mount's position log has
  > it standing still at altitude 70 deg

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

Run over the archive on 2026-10-20 (headers and solutions) and 2026-10-21 (the position logs), read-only. Every
episode is listed in `docs/telescope_faults_census.csv` (306 rows: telescope, night, kind, target, UT range, frames,
minutes, detail). Scripts and raw output are in `mh_scratch/telescope_faults/` (`census_tracking.py`,
`census_drift.py`), `mh_scratch/sitelog_tracking/` (`census_sitelog.py`, `census_fix.py`, `crosscheck_track.py`,
`trail_sample.py`, `*.jsonl`) and, for Artemis, `mh_scratch/artemis_tracking/` (`acp_frames.py`,
`census_artemis.py`, `clock_check.py`, `final_offsets.py`, `verify_frames.py`, `crosscheck_artemis.py`,
`trail_artemis.py`, `census_final.jsonl`) on appct.

### Tracking stops

The census read the Astra-era raw headers: 1,615 telescope-nights with light frames and 2.77 million frames. The SSO
raw names carry each exposure's start, so frames were sampled every 5 minutes without being opened. Each run of stopped
samples then had its first and last stopped frame found by bisection: 171,000 header reads, where reading every
header beside the queue's jobs on the busy raw-data array would have taken days. A stop shorter than 5 minutes can
fall between samples. One false stop, on Callisto 2024-12-09, is not counted: TRACKING was 0 for half an hour, but the airmass
moved as only a tracking mount's does.

| Year | Io | Europa | Ganymede | Callisto |
|---|---|---|---|---|
| 2023 | | | | 1 / 11 |
| 2024 | | | | 6 / 295 (4 h) |
| 2025 | | | 7 / 242 (9 h) | 1 / 275 |
| 2026 | 0 / 92 | 0 / 215 | 2 / 237 (3 h) | 0 / 248 |

*Telescope-nights with a stop / Astra-era nights read (hours of frames lost).*

Seventeen stops on 17 nights, of three kinds:

- **Six long stops mid-sky**, 1,353 frames and about 14 hours in all:
  - Ganymede 2025-07-16 TOI-4336 c, 570 frames
  - Ganymede 2025-08-25 Sp1954+1801, 333 frames
  - Ganymede 2025-09-04 Sp1954+1801 again, 02:59–04:38 UT, 133 frames
  - Ganymede 2026-05-05 Sp1649+0444, 21 frames
  - Ganymede 2026-06-30 Sp2253-4801, 103 frames
  - Callisto 2024-10-12 Sp0002+0115, 02:04–04:22 UT, 193 frames
- **Six at the mount's altitude limit.** The mount stopped at 20 degrees (airmass 2.898) while the camera kept
  exposing. These were Chilean and asteroid targets, and Sp0254-1934 on Callisto twice (2024-11-19 and 11-27, 12 and
  27 frames).
- **Five single frames**, with the mount idle for 52–85% of the exposure. The position log has the mount tracking
  through all five: each is the first frame of a block, and the hour-angle axis already turned at the sidereal rate
  (Callisto 2023-12-28: the slew ended 15 s before the exposure began). Astra's flag lags at block starts.

The position log, read the same way over the Astra era (ESO's frames, 1,594 nights the log covers), finds 18 stops.
They are the twelve stops above other than the five single frames, with one more frame on four of them, and six the
header census lacks. Two are on disk and shorter than its 5-minute sampling, and Astra's `TRACKING` agrees on both (Callisto
2024-02-29 ch_TOI-519b, Ganymede 2025-08-22 ch_HS2231+2441_2, both at the 20-degree limit). Four are on Chilean
programmes' frames, which are never downloaded (Io 2026-02-07 and 06-03/04, Ganymede 2025-07-21, at 20–22 degrees).
On Callisto 2024-12-09, the false stop above, the log has the mount tracking.

### Tracking stops before Astra

The position log against ESO's list of every SCIENCE frame (the completeness ledger's per-frame inventory,
`mh_scratch/ledger/eso/frames`, 2026-10-08: start, exposure, object, programme), so without opening a frame:
7,002 MaxIm DL telescope-nights from 2019 to Astra, 5.5 million frames, the log covering 5.2 million of them on 6,695
nights. Each log was read only between the night's first and last frame, a night a second at most. Every night where
some frame overlapped a still period was read again through `MountLog` itself, with the clock set from ACP's slews,
and those results stand.

| Year | Io | Europa | Ganymede | Callisto |
|---|---|---|---|---|
| 2019 | 1 / 259 (6 h) | 1 / 271 (5 h) | 1 / 263 (6 h) | 1 / 284 (1 h) |
| 2020 | 0 / 155 | 0 / 158 | 0 / 146 | 0 / 143 |
| 2021 | 0 / 278 | 0 / 324 | 0 / 319 | 25 / 293 (68 h) |
| 2022 | 0 / 287 | 4 / 286 (6 h) | 8 / 257 (21 h) | 7 / 221 (21 h) |
| 2023 | 2 / 316 | 3 / 318 (6 h) | 3 / 302 | 0 / 251 |
| 2024 | 0 / 289 | 0 / 271 | 1 / 299 (1 h) | |
| 2025 | 0 / 294 | 3 / 323 (4 h) | 1 / 47 (1 h) | |
| 2026 | 0 / 19 | 0 / 22 | | |

*Telescope-nights with a stop / MaxIm DL nights the log covers (hours of frames lost).*

79 stops on 61 telescope-nights, 15,331 frames and 147 hours; 74 of them, 143 hours, on the survey's own programme.
Fifty stops of an hour or more hold 141 of the hours, and 41 lasted to the night's last frame. Most came in spells:

- **Callisto, May 2021 to November 2022.** 37 stops on 32 nights, 89 hours. Through 2021 and into March 2022 the
  controller logged errors on axis 0 (`ERR_DeviceError`), mostly keeping a trajectory running while the axes stood
  still (2021-05-31 Sp2038-5223, 230 frames to the end of the night; 2021-06-09 and 06-10 Sp1609-3431, 367 and 531;
  2021-11-18 Sp0335-4430, 675). From May 2022 its GPS clock lost sync (`GPS clock is not synchronized`, `Time not
  valid. tracking points deleted`) on each night the mount stopped: 05-11 Sp1359-5636 (982 frames), 05-19
  Sp1424-3514 (1,085), 10-24 Sp0027-5401 (746), 11-09 Sp0218-0617 (625), 11-14 Sp0339-3525 (1,046).
- **Ganymede, 2022-02-16 to 03-01.** 16 stops on 8 nights, 22 hours, with `ERR_DeviceError` on axis 0 and, on 03-01,
  129 `axis (0) unexpectedly changed to powered on state`. Sp0830+0947 lost 437 and 343 frames on 02-18 and 02-19.
- **2019-05-05, all four telescopes.** Io, Europa and Ganymede stopped mid-target at 00:23–00:44 UT and again at
  05:24–06:10, Callisto once at 00:50; Io and Ganymede logged `ERR_TrajectoryNoPoints` (`Trajectory ran out of
  points`). The second stops ended at 07:14–07:21 UT, as the controllers were shut down (`TSI is closing down`): 7
  stops, 1,739 frames, 17 hours.
- **Europa, 2022 to 2025.** Ten nights, one stop each, mostly from 06:15–09:30 UT to the end of the night (2023-09-20
  Sp0216+1335, 548 frames; 2025-04-30 Sp1621-3711, 206).

Three stops were at the altitude limit, here 23.5–24 degrees (Io 2023-05-02 and 06-26 TOI-1227.01, 11 frames each;
Europa 2022-02-01, a Chilean target).
Not counted: Europa 2019-09-04 Gaijin-1 (157 frames), a manual observation after ACP's script was stopped at
06:26 UT, with a trajectory running throughout, so probably deliberate; and four 1-s TEST exposures on engineering
nights. SSO nights before February 2019 have no position log. Artemis's frames are not at ESO, so it has a census of
its own (next section).

### Tracking stops on Artemis

Artemis (Teide, MaxIm DL throughout) has position logs from 2021-08-25 (regularly from 09-01) to 2025-06-10, but
they hold only their header from 2022-06-24 to 2023-07-28. Its frames are not at ESO, so there is no inventory like
SSO's, and its raw names (`Sp1356+4342-S001-R001-C012-I+z.fts`) carry no time. ACP's log, beside the position log,
has what the census needs. Each plan image is an `Imaging to <name>` line, the raw frame's name, followed by
`(taking N sec. exposure; <filter> filter; ...)` and `(starting exposure)`, by the ACP PC's clock, which also stamps
`DATE-OBS`. Pointing, focus and calibration exposures have no `Imaging to` line. Checked against the raw headers on 7
nights from 2021 to 2025 (182 frames read), `DATE-OBS` fell within 0.65 s of ACP's second, and `EXPTIME` and `OBJECT`
agreed on every frame. Every light frame on disk had its ACP lines on those nights. When a target's plan is
restarted, ACP reuses the image names, and the later frames overwrite the earlier ones (Cosmos_2539 on 2025-03-20).

The census read the 936 nights with a position log. Artemis took 482,153 light frames on 646 of them, and the log
covers 475,959. ACP's log lists 481,228 of the frames; on four nights its copy misses part of the night, and the 925
frames on disk that it lacks were read from their headers instead. The log covers nothing on the 292 nights with
frames in the header-only spell (218,000 frames) or on 49 nights with no position log at all (37,000 frames, 12 of
the nights in October 2024). Artemis's nights before 2021-08-25 and after 2025-06-10 have no position log either, so
T11 says `unknown` for them unless T7 finds trailed frames. On 9 of the covered nights (2024-08-02 to 2025-02-05) ACP
took 4,609 frames that are not on disk: those nights have no images directory.

**The clock.** Artemis's position log runs 1–3 s ahead of ACP's clock, with no spell like Europa's in 2019.
`clock_offset()` finds 1.0–3.0 s on 396 of the 646 nights (median 2.0 s). On 212 it has too few slews and returns
nothing, so T11 reads those nights with the offset 0. On 38 nights (2021–2025, 17 of them in 2024) it returns −42 to
−45 s, which is wrong. Artemis focuses at the start of each block. FocusMax slews to a focus star and back, which
the position log shows and ACP's log has no `Start slew to` line for. ACP then takes a pointing exposure and slews
to the target again, a move of about an arcminute that the log does not show. Its `Start slew to` line therefore
comes 43 s after FocusMax's return slew. This happens 4–5 minutes after the acquisition, so the 150-s exclusion of
a second slew to the same target does not catch it, and on nights with two or three targets these pairs outvote
the acquisitions. The census pairs instead each ACP slew that took 15 s or more, from `Start slew to` to `(slew
complete)`, with the longest motion the log shows within 150 s of it, and takes the median of the pairs within
−5 to +10 s, a far pair being a mispairing. That gives +1 to +3 s on all 38 nights, and +1 to +5 s on 151 of the
212 nights without an offset. The census uses those offsets, and 0 on the other 61. Read with
`clock_offset()`'s own offsets, as T11 reads them now, every night has the same frames stopped, so the wrong
offsets change no report here.

| Year | Artemis |
|---|---|
| 2021 | 1 / 90 (1 h) |
| 2022 | 0 / 107 |
| 2023 | 2 / 74 (4 h) |
| 2024 | 0 / 242 |
| 2025 | 1 / 133 |

*Nights with a stop / Artemis nights the log covers (hours of frames lost).*

Four stops on four nights, 691 frames on disk and 5 hours, all at one limit of the mount's hour-angle axis. In the
646 nights the log shows no other stop, nothing like the drive failures on Callisto and Ganymede. Four single
frames, none of them saved, are not faults.

- **The hour-angle limit.** All four stops were on targets far north: TOI-2267.02 and .03 at declination +85, and
  TOI-2096 at +75. Each block began with the target on the far side of the pole, 6.9–8.3 h east of the meridian,
  and the mount's hour-angle axis past 180 degrees (`HA(inst.)` 236–257 in the log). The axis turned at the
  sidereal rate until it read 265.0 degrees, the target then 6.3 h east, and stopped there, 0.5–2 h after the
  acquisition, while ACP kept exposing:
  - 2021-09-12 TOI-2267.02, 23:39–00:45 UT, 197 frames to the block's end
  - 2023-10-12 TOI-2267.02, 21:45–00:16 UT, 449 frames to the block's end
  - 2023-10-31 TOI-2267.03, 20:31–21:51 UT, 38 frames to the block's end
  - 2025-01-14 TOI-2096 (`_TOI2096` in the frames), 21:14–21:17 UT, 8 frames, until ACP's script was interrupted
    at 21:18 UT. The 8th frame was cut short and not saved. The next run acquired the target again at 21:23 with the
    axis on the other side of its range (−65 degrees at the block's end), and the mount tracked to the end.

  On 2021-09-12 and 2025-01-14 the controller logged `ERR_TrajectoryNoPoints` (`Trajectory ran out of points`) as
  the mount stopped; the 2023 stops left no event. ACP's plans do not know the limit. At declination +85 a stopped
  mount trails stars by only 1.3 arcsec a second, so the 10-s frames of TOI-2267.02 show streaks of about 40 pixels.
- **Not faults.** Three single frames at a restart of ACP: 2022-04-09 Sp1256-1257, 2024-12-03 Sp0320+1854 and
  2025-01-29 TOI-7185.01. ACP's script was interrupted during the exposure (`The script was interrupted before
  completion`), and tracking stopped until the next run turned it on again; on 2024-12-03 and 2025-01-29 the plan
  that followed turned it off (`RUN STOPPED AT SPECIFIED QUIT TIME`, `turning tracking off for safety`). MaxIm did not
  save these frames, and on 2024-12-03 and 2025-01-29 the next run gave the same name to a later frame. The fourth
  was the last frame of 2023-10-08: 56 s into its 180 s the controller ended its trajectory (`Unexpected end of
  trajectory execution`, `Tracking position bad`) and the mount went to its park position as the weather turned,
  which ACP then reported (`Weather is unsafe!`). That frame was not saved either.

Checked as for SSO. `CURRENT.TRACK` in the controller's status log was 0, no trajectory running, from the first
frame to the last of all four stops. During the single frames the trajectory ended (2023-10-08, 2024-12-03) or none
ran (2022-04-09), and on 2025-01-29 one started 148 s in, as ACP turned tracking on again. By their `DATE-OBS`,
691 of the four stops' 692 frames are on disk as the exposures ACP logged (`verify_frames.py`); the other is
2025-01-14's 8th. In the frames
(`trail_artemis.py`: the frame before each stop, its first, middle and last frames, and the frame after), the stars'
axis ratios are 2.9–17 inside the stops and 1.0–1.2 outside. The exception is the frame before the 2025-01-14 stop,
at 2.5, during which the mount stopped. T11's own path, `night()` over the raw headers, finds the same 197 frames on
2021-09-12.

### Pointing drifts

The census read the per-target `*_output.fits` imagelists: each photometered frame's WCS, one read per target-night,
from v2 and, where v2 has none, from v3. That covered 18,177 target-nights on 9,908 telescope-nights.

| Year | Io | Europa | Ganymede | Callisto | Artemis |
|---|---|---|---|---|---|
| 2017 | 14 / 46 (10 h) | 16 / 192 (4 h) | | | |
| 2018 | 6 / 182 (5 h) | 8 / 279 (29 h) | 0 / 55 | 12 / 167 (9 h) | |
| 2019 | 0 / 280 | 0 / 284 | 3 / 295 (2 h) | 0 / 308 | 0 / 162 |
| 2020 | 0 / 148 | 0 / 154 | 0 / 141 | 0 / 139 | 1 / 213 |
| 2021 | 0 / 269 | 0 / 312 | 0 / 292 | 0 / 273 | 2 / 266 (1 h) |
| 2022 | 0 / 256 | 1 / 270 | 2 / 255 (6 h) | 3 / 56 (14 h) | 1 / 228 (2 h) |
| 2023 | 2 / 265 | 0 / 274 | 0 / 269 | 0 / 195 | 0 / 181 |
| 2024 | 1 / 246 (3 h) | 0 / 222 | 0 / 265 | 8 / 194 (1 h) | 6 / 154 (9 h) |
| 2025 | 0 / 229 | 7 / 266 (11 h) | 12 / 254 (4 h) | 0 / 212 | 1 / 201 (3 h) |
| 2026 | 0 / 86 | 1 / 217 (1 h) | 0 / 212 | 5 / 226 (3 h) | 0 / 137 |

*Telescope-nights with a drift / telescope-nights with a measured pointing (hours affected).*

There are 195 episodes on 112 telescope-nights. In 71 of them, on 51 nights, the target left the field, and 108
arrived as jumps.

- **2017–18 commissioning** (Io, Europa, Callisto; 56 nights). The telescope pointed up to 3 degrees away for minutes
  to hours, with the frames still named for the target and the target off the field.
- **Artemis, February–March 2024** (6 nights). There were repeated jumps of 2.6–7.6 arcmin, with the target off the
  field on four of the nights: TOI-5720 on 02-10 and 02-15, Sp1356+4342 on 02-11, Sp1444+3002 on 02-14, and
  Sp1438+6408 on 03-10 and 03-13.
- **Ganymede and Europa in 2025** (19 nights, 21 episodes). Fifteen start with their block and last 7–48 minutes at
  1.2–3.5 arcmin: the telescope settling for longer than the 10 minutes allowed. The others are mid-block, such as
  Europa 2025-03-20 TOI-1227, 2.8 arcmin over five hours from 04:54 UT.
- **Callisto 2026-10-14 to 10-16, Sp0439+1615.** Each night from 04:30 to 05:30 UT the target sat just off SPIRIT's
  field. The pointing then moved 1.7 arcmin and put it on the field.

The census sees photometered frames only. Frames that T7 or the photometry lost are not measured: on Artemis
2024-02-11, nothing after 01:49. A rerun with the known-scale solver (v4.1.0) measures more of them.

221 target-nights could not be judged, because their stored solutions scatter from frame to frame (163 of them are
Callisto 2024). Those products were made with wrong WCSs on many frames. This is not a telescope fault but a
candidate list for reruns: `mh_scratch/telescope_faults/noisy_wcs_target_nights.csv`.

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
| `STILL_RATE` | 0.1 | share of the sidereal rate below which the position log has an axis standing still |
| `LOG_STEP_MAX` | 10 s | lines of the position log further apart leave the time between them unknown |
| `LOG_COVER_MIN` | 0.5 | share of an exposure the position log must cover to say whether the mount tracked |
| `SLEW_RATE` | 0.05 deg/s | an axis faster than this in the position log is slewing |
| `OFFSET_MATCH`, `OFFSET_AGREE`, `OFFSET_VOTES` | 150 s, 3 s, 2 | ACP's slews and the log's within 150 s pair; 2 or more agreeing to 3 s set the clock offset |
