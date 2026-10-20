# Shutter detection

`calibration/shutter_detect.py` measures, for every night, whether the SSO Andor cameras were read out with or
without the mechanical shutter, and at which readout rate. The pipeline has to know this itself: CAMSHUTT is a
static value in Astra's config, changed by hand, so it would not change if a shutter failed open, and
`shutter_overrides.csv` only covers the windows somebody has verified.

## Order of evidence

First that gives an answer wins (agreed 2026-10-19):

1. **Override table** (`calibration/shutter_overrides.csv`): human-verified windows. Applied frame by frame in
   `smear.shutter_state()`; the night's record notes any window the measurement disagrees with.
2. **The night's own flats**: the frames physically in `Observations/<TEL>/images/<NIGHT>/`, never the reduction's
   flat lists (which borrow from earlier nights: Ganymede 20250408-09 was reduced with 20250407's shuttered morning
   flats while it was shutterless). Each flat gets the readout-ramp test. Flats before the first science frame are the
   evening twilight, after the last the morning. When both twilights agree, the state holds for the whole night.
3. **Science frames** of the night, where the flats are missing, one-sided or disagree: the star-trail excess,
   stacked over 24 frames of an exposure time up to 40 s, plus bins of 12 frames at both ends of the night (a
   whole-night stack would average a switch away). On a night whose ends differ, bins bisect the night to bracket the
   switch. For longer exposures, on a night nothing else answered, the sky readout ramp against a shuttered reference
   flat (from the nearest earlier night within 10 days) is measured too; only a 1 MHz reading of it is used, and only
   as probable (see Validation).
4. **CAMSHUTT**, when the measurement is inconclusive.
5. **Shutter used.**

The record is `reduction/shutter.json`, written at T2 (pipebias), read back at T4 and T6. A night borrowed for its
flats gets `shutter_<night>.json` at T4, measured on that night's flats only. `python -m calibration.shutter_detect
<night dir> [--all-levels] [--json]` measures a night and writes nothing.

## What the measurement does

Every calibrated frame (T4 flats, T6 science) carries:

| card | meaning |
|---|---|
| `SHUTDET` | measured state: `SHUTTERED`, `SHUTTERLESS` or `UNKNOWN` |
| `SHUTDSRC` | evidence: `flats`, `flats+science`, `science`, `sky_ramp` |
| `SHUTDCNF` | `unambiguous`, `probable` or `inconclusive` |
| `SHUTDTRD` | measured readout time (shutterless only) |
| `SHUTDIS` | T when the measurement disagrees with the frame's CAMSHUTT |
| `SMEARSRC` | what decided the de-smear (now also on shuttered frames that have a record) |

Disagreements with CAMSHUTT, an override window or READOUTM, and shutterless nights that no header or window
describes, are recorded in `shutter.json`, printed as warnings in the T2/T4/T6 logs, and listed in the notifier's
morning digest ("Shutter state measured against the headers", nights reduced in the last 24 h).

**The measurement changes what is de-smeared only where `"detect_auto_apply": true` is set under `"smear"` in
`instrument_config.json`, and only for unambiguous segments.** It is `true` for all four SSO telescopes since
2026-10-20 (decisions below). Where it is off, the decision stays override, CAMSHUTT, default, and the measurement
is reported only.

## The rules

Each rule below makes a segment **unambiguous**, i.e. eligible for auto-application. Everything else is
`probable` or `inconclusive` and is reported only. The numbers are the module constants in `shutter_detect.py`.

| # | rule | evidence behind it |
|---|---|---|
| A | **Each flat**: shutterless if its ramp t_read > 0.75 s and the ramp it implies > 5% across the frame; shuttered if \|t_read\| < 0.5 s and a 1.4 s readout would show > 3% (flats up to ~33 s); else inconclusive, as are flats below 2,000 ADU or above 85% of the mode's saturation | `smear_check.flat_verdict`; 1,653 SSO flats (2025-26) separated cleanly |
| B | **A twilight** is conclusive with ≥ 2 conclusive flats, all agreeing. Its readout rate is the intercept of t_read against exposure time (Theil–Sen), within 15% of a measured rate (1 MHz 4.33 s, 3 MHz 1.5 s); a twilight straying > 25% + 0.1 s from its line mixes rates | apparent t_read grows ~0.009 s per s of exposure on shuttered flats (flat-field structure) and the twilight adds or subtracts; intercepts: shuttered 0.00–0.02 s, 1 MHz 4.29 s evening / 4.42–4.46 s morning, 3 MHz 1.47–1.59 s |
| C | **Both twilights agree** (state and rate): the whole night | the agreed hierarchy |
| D | **One twilight confirmed by the science**: a shutterless twilight by a conclusive trail at its rate; a shuttered one when every trail stack (the whole-night samples and the end-of-night bins) rules out 1 MHz | a shuttered stack is not required to exclude 3 MHz: the trail test cannot, a short shuttered flat already has |
| E | **A split located by the science between two disagreeing twilights**: both sides; the bracket between them stays with CAMSHUTT | bisection with 12-frame bins, never puts the switch outside its bracket (unit test over all positions) |
| F | **Science alone** (no usable flats): **never, by default** (`SCIENCE_ALONE_SIGMA = inf`). A finite value makes a conclusive trail with that margin (σ) to every other hypothesis unambiguous | the trail bias has a long tail: Ganymede 20230130, shuttered, read −5.11 ± 0.51 s at 31 s; a positive bias that size would read as a conclusive 1 MHz and de-smear a whole shuttered night |
| — | The **sky ramp** gives at most a probable answer, and only when it reads 1 MHz; a **READOUTM** that disagrees with the measured rate makes the segment probable | shuttered nights read +1.0 to +2.6 s on it (3 MHz) |

Star trails (rules D–F) are judged with a PSF-asymmetry allowance of max(0.5 s, 0.02 s × exposure) added in
quadrature to the statistical error. Shutterless rates are intervals, 0.9–1.35 × their time, because a broad PSF
inflates the trail. Hypotheses are consistent within 3σ and excluded beyond it. Without READOUTM the hypotheses are
shuttered, 1 MHz and 3 MHz; with READOUTM on every frame, shuttered and the rate it names.

**Decisions (Matthew, 2026-10-20), in force from v4.2.0:**

1. **Rules C, D and E are applied** (two twilights agree; one twilight confirmed by the science; a split between
   disagreeing twilights). These carried every unambiguous answer in the validation, all correct.
2. **Rule F is not**: a night is never settled by its science alone (`SCIENCE_ALONE_SIGMA = inf`). The three
   verified science-only shutterless nights (Ganymede 20250408, 20250409, 20260216; margins 4.0–4.9σ) are in the
   override table; the −5.1 s bias above says a shuttered night could one day pass. Such a night costs ~1% noise until
   someone adds an override window, which the digest prompts.
3. **All four SSO telescopes** (`"detect_auto_apply": true` under each `"smear"`). Artemis has no smear model.
4. **3 MHz without flats is accepted as it is.** Without READOUTM the trail test cannot tell 3 MHz from shuttered
   (both sit within the PSF bias), so a 3 MHz night without both twilights stays probable or inconclusive and falls
   back to CAMSHUTT. Astra writes READOUTM from 2026-10-18, which narrows the hypotheses for new frames.

## Validation

Run 2026-10-19 on appct (orchard-server, read-only on production; scripts and records in
`mh_scratch/shutter_detect/` on appcs).

**Sample.** 116 night directories: all 34 Ganymede nights from 2025-04-04 to 04-12 and 2026-02-10 to 03-06 (the three
verified shutterless windows and their neighbours), and on each SSO telescope 15 random nights from 2025-26 and 5 from
before 2025. Eight directories in the sample hold no frames, leaving **106 nights**: 17 in the verified shutterless
windows, 71 shuttered nights from 2025-26 (verified: the 1,653-flat survey found no other shutterless period) and 18
presumed shuttered from before 2025. Truth per night is the override window holding <night>+1 03:00 UTC (every
verified switch happened in daytime). Every test was run on every night (`--all-levels`) to score each level on its
own; the hierarchy's answers come from a second run in production mode with the final code.

**Headline.** No night was answered wrongly at any confidence, and no flat, twilight, trail stack or bin on a
shutterless night was read as shuttered.

| truth | nights | unambiguous | probable | inconclusive | wrong |
|---|---|---|---|---|---|
| shutterless (verified windows) | 17 | 9 | 6 | 2 | 0 |
| shuttered, 2025–26 (verified) | 71 | 45 | 18 | 8 | 0 |
| shuttered, before 2025 (presumed) | 18 | 11 | 6 | 1 | 0 |
| **all** | 106 | 65 | 30 | 11 | 0 |

Unambiguous answers came from both twilights agreeing (51) or one twilight confirmed by the science (14); none
needed the science alone. The inconclusive nights have no usable flats of their own (none, or one) and science the
trail test cannot judge (short exposures with too few bright isolated stars, or shuttered trails that cannot exclude
3 MHz), or no science at all (Ganymede 20260215 and 20260217, calibrations only, inside the February window). They
fall back to CAMSHUTT or the default, which is right for the shuttered ones.

**What the levels can and cannot do.**
- *Flats* separate cleanly when judged per twilight: per flat, 8 of 1,965 shuttered flats read as shutterless, all
  among the first, longest flats of two dark early-morning twilights (Callisto 20250218, Ganymede 20260308: 3–28 s
  at 27–60 s, falling to ~0 by 20 s). Both twilights came out *mixed*, so neither was used. A twilight's t_read rises
  linearly with exposure: the intercept is the readout time (shuttered 0.00–0.02 s, 1 MHz 4.29 s evening and
  4.42–4.46 s morning, 3 MHz 1.47–1.59 s), whereas per flat 1 MHz morning flats reach 5.5 s at 60 s.
- *Star trails* see 1 MHz at every exposure up to 40 s (8 of 8 shutterless stacks, margins 3.7–6.0σ), but not
  shuttered from 3 MHz: the PSF bias on shuttered nights spans −5.1 to +2.8 s (+7.2 ± 2.9 s once). So a trail confirms
  a shutterless twilight, rules out 1 MHz for a shuttered one, and locates a switch; it never sets a night on its own.
- *Sky ramp* (exposures > 40 s): not reliable. Shuttered nights read +1.0 to +2.6 s against a shuttered reference
  (and −0.9 s once): night sky and twilight light the detector differently by ~1%. Four of 69 shuttered stacks were
  judged 3 MHz; the 1 MHz nights read 5.2–5.4 s at 41 s. Only a 1 MHz reading is used, as probable.

**End to end.** T2, T3, T4 and T6 (8 frames) were run in scratch with the live release (20261019T121959Z-e585ad8)
and with this branch on Ganymede 20260221 (shutterless, override window) and Europa 20261012 (shuttered,
CAMSHUTT T, borrowing 20261011's flats). All 11 and 12 FITS products (master bias, darks, flats, calibrated frames)
have identical pixels; the only header differences are the SHUTD* cards and SMEARSRC on shuttered frames. Each run
wrote `shutter.json`, and Europa's T4 wrote `shutter_20261011.json` for the borrowed night.

**Headers.** No disagreement was recorded on any night: the measurement agreed with all three override windows,
and with CAMSHUTT on the three sampled nights that carry it (from 2026-07-17). No sampled frame carries READOUTM
(Astra writes it from 2026-10-18), so the READOUTM cross-check, the CAMSHUTT disagreement path and the digest lines
are covered by the unit tests, and by Europa 20261012 in the end-to-end run (CAMSHUTT T, SHUTDIS F).

### Each level on its own (`--all-levels`)

FP = read as shutterless on a shuttered night; FN = read as shuttered on a shutterless night; wrong rate = shutterless at the wrong readout rate.

| level | truth | n | correct | FP | FN | wrong rate | inconclusive |
|---|---|---|---|---|---|---|---|
| each flat (ramp test) | shuttered | 1965 | 1635 (83%) | 8 (0.4%) | 0 | 0 | 322 (16%) |
| each flat (ramp test) | shutterless | 319 | 271 (85%) | 0 | 0 | 0 | 48 (15%) |
| each twilight (≥ 2 flats) | shuttered | 132 | 123 (93%) | 0 | 0 | 0 | 9 (7%) |
| each twilight (≥ 2 flats) | shutterless | 18 | 18 (100%) | 0 | 0 | 0 | 0 |
| trail, whole-night stack (24 frames) | shuttered | 99 | 13 (13%) | 0 | 0 | 0 | 86 (87%) |
| trail, whole-night stack (24 frames) | shutterless | 8 | 8 (100%) | 0 | 0 | 0 | 0 |
| trail, end-of-night bin (12 frames) | shuttered | 128 | 18 (14%) | 0 | 0 | 0 | 110 (86%) |
| trail, end-of-night bin (12 frames) | shutterless | 14 | 11 (79%) | 0 | 0 | 0 | 3 (21%) |
| sky ramp (experimental) | shuttered | 69 | 6 (9%) | 4 (6%) | 0 | 0 | 59 (86%) |
| sky ramp (experimental) | shutterless | 9 | 3 (33%) | 0 | 0 | 0 | 6 (67%) |

### The hierarchy's answer per night, by confidence

From the production-mode run (final code, no `--all-levels`), 106 of 106 nights.

| truth | confidence | n | correct | FP | FN | wrong rate | inconclusive | split (wrong) |
|---|---|---|---|---|---|---|---|---|
| shuttered | unambiguous | 56 | 56 | 0 | 0 | 0 | 0 | 0 |
| shuttered | probable | 24 | 24 | 0 | 0 | 0 | 0 | 0 |
| shuttered | inconclusive | 9 | 0 | 0 | 0 | 0 | 9 | 0 |
| shutterless | unambiguous | 9 | 9 | 0 | 0 | 0 | 0 | 0 |
| shutterless | probable | 6 | 6 | 0 | 0 | 0 | 0 | 0 |
| shutterless | inconclusive | 2 | 0 | 0 | 0 | 0 | 2 | 0 |

Trail PSF bias on shuttered nights (whole-night stacks), apparent t_read in s: Callisto -0.67..+0.35 (n=9); Europa -1.23..+2.79 (n=19); Ganymede -5.11..+7.23 (n=26); Io -1.97..+2.83 (n=19).

No shuttered trail stack or bin was ever conclusive at 1 MHz. The nearest was a 12-frame bin on Europa 20250331 (30 s) at 5.43 ± 1.96 s: inside the 1 MHz interval, but within 2.7σ of zero, so inconclusive.
Margins of the conclusive whole-night 1 MHz stacks on shutterless nights: 3.7 sigma (Ganymede_20260220 39 s), 3.8 sigma (Ganymede_20260219 39 s), 3.9 sigma (Ganymede_20260218 39 s), 4.0 sigma (Ganymede_20260221 39 s), 4.0 sigma (Ganymede_20260216 39 s), 4.5 sigma (Ganymede_20250409 20 s), 4.9 sigma (Ganymede_20250408 20 s), 6.0 sigma (Ganymede_20250409 10 s).

Production-mode run time per night (header scan, flats, and the science only where the flats do not settle the night): median 22 s, 90th percentile 60 s, max 182 s, over 106 nights.

Answers that differ between the all-levels run and the production-mode run: none.

### Synthetic split nights

Real frames, linked into one night directory: Ganymede 20250407 (shuttered) followed by 20250409 (shutterless, 1 MHz), 20 s exposures; the true switch lies between 2025-04-08T07:15:56 and 2025-04-10T03:58:43.

| case | result | bracket found | brackets the true switch | segments |
|---|---|---|---|---|
| 20250407's own flats first | split (flats+science, probable) | 2025-04-08T07:07:21 to 2025-04-10T04:05:28 | yes | 04-07T23:11–04-08T07:07 shuttered (unambiguous); 04-08T07:07–04-10T04:05 unknown (inconclusive); 04-10T04:05–04-10T07:17 shutterless 1MHz (probable) |
| no flats | split (science, probable) | 2025-04-08T07:07:21 to 2025-04-10T04:05:28 | yes | 04-07T23:11–04-08T07:07 unknown (inconclusive); 04-08T07:07–04-10T04:05 unknown (inconclusive); 04-10T04:05–04-10T07:17 shutterless 1MHz (probable) |

### Every night

Flats: conclusive/all, state, rate (intercept t_read). Science: whole-night trail stacks and sky ramps (t_read in s), end bins omitted. Truth from the verified override windows.

| telescope | night | truth | evening flats | morning flats | science | result | source | confidence | disagreements |
|---|---|---|---|---|---|---|---|---|---|
| Callisto | 20180709 | shuttered | 10/10 shuttered | 10/10 shuttered | trail 25s 0.22 | shuttered | flats | unambiguous | – |
| Callisto | 20211129 | shuttered | – | 10/10 shuttered | trail 13s -0.09 | shuttered | flats+science | unambiguous | – |
| Callisto | 20220127 | shuttered | 7/7 shuttered | 7/7 shuttered | trail 12s -0.37 sh; sky 69s -0.22 sh | shuttered | flats | unambiguous | – |
| Callisto | 20240825 | shuttered | 20/20 shuttered | 15/17 shuttered | trail 17s 0.35; trail 40s -0.67 | shuttered | flats | unambiguous | – |
| Callisto | 20250127 | shuttered | 20/20 shuttered | 7/13 shuttered | trail 15s 0.23; trail 16s -0.01 | shuttered | flats | unambiguous | – |
| Callisto | 20250218 | shuttered | 20/20 shuttered | 8/15 mixed | trail 10s -0.24 sh; trail 14s -0.25 sh | shuttered | flats+science | unambiguous | – |
| Callisto | 20250426 | shuttered | – | 9/15 shuttered | trail 7s –; trail 8s – | shuttered | flats | probable | – |
| Callisto | 20250728 | shuttered | 15/15 shuttered | 11/15 shuttered | trail 3s –; trail 8s – | shuttered | flats | unambiguous | – |
| Callisto | 20251014 | shuttered | 20/20 shuttered | 15/20 shuttered | trail 9s –; trail 10s – | shuttered | flats | unambiguous | – |
| Callisto | 20251018 | shuttered | 0/6 inc. | 15/20 shuttered | trail 10s –; trail 25s – | shuttered | flats | probable | – |
| Callisto | 20251029 | shuttered | 20/20 shuttered | 15/20 shuttered | trail 2s –; trail 25s – | shuttered | flats | unambiguous | – |
| Callisto | 20251113 | shuttered | 32/35 shuttered | 16/26 shuttered | trail 14s –; sky 45s – | shuttered | flats | unambiguous | – |
| Callisto | 20251231 | shuttered | 30/30 shuttered | 23/30 shuttered | trail 6s –; trail 12s – | shuttered | flats | unambiguous | – |
| Callisto | 20260313 | shuttered | 30/30 shuttered | 23/30 shuttered | trail 6s –; trail 11s – | shuttered | flats | unambiguous | – |
| Callisto | 20260331 | shuttered | 25/30 shuttered | 0/5 inc. | trail 13s –; trail 16s – | shuttered | flats | probable | – |
| Callisto | 20260803 | shuttered | – | – | trail 6s –; trail 12s – | unknown | – | inconclusive | – |
| Callisto | 20260919 | shuttered | 9/26 shuttered | 11/17 shuttered | trail 9s –; trail 17s – | shuttered | flats | unambiguous | – |
| Callisto | 20261007 | shuttered | – | 12/17 shuttered | trail 6s –; trail 15s – | shuttered | flats | probable | – |
| Callisto | 20261010 | shuttered | 17/30 shuttered | – | trail 6s –; trail 15s – | shuttered | flats | probable | – |
| Europa | 20170318 | shuttered | – | – | sky 120s –; sky 120s –; sky 120s – | unknown | – | inconclusive | – |
| Europa | 20171102 | shuttered | 6/6 shuttered | 10/10 shuttered | trail 12s 0.56; trail 23s -0.18 | shuttered | flats | unambiguous | – |
| Europa | 20220122 | shuttered | – | 7/7 shuttered | trail 10s 0.37; trail 16s -0.38 sh | shuttered | flats+science | unambiguous | – |
| Europa | 20220519 | shuttered | 7/7 shuttered | 5/7 shuttered | trail 13s 0.18; trail 26s 0.21 | shuttered | flats | unambiguous | – |
| Europa | 20221112 | shuttered | 12/12 shuttered | 12/12 shuttered | trail 10s -0.33 sh; trail 20s -1.23 sh; sky 63s 0.14 | shuttered | flats | unambiguous | – |
| Europa | 20250125 | shuttered | – | 12/12 shuttered | trail 26s -0.81 sh; trail 39s – | shuttered | flats+science | unambiguous | – |
| Europa | 20250331 | shuttered | – | 12/12 shuttered | trail 30s 0.91 | shuttered | flats | probable | – |
| Europa | 20250423 | shuttered | 12/12 shuttered | 12/12 shuttered | trail 17s 0.06; trail 21s 0.44 | shuttered | flats | unambiguous | – |
| Europa | 20250805 | shuttered | 7/7 shuttered | 7/7 shuttered | sky 82s 0.12; sky 97s 0.42 | shuttered | flats | unambiguous | – |
| Europa | 20250807 | shuttered | 7/7 shuttered | 7/7 shuttered | sky 82s 0.82; sky 97s 2.40 | shuttered | flats | unambiguous | – |
| Europa | 20250920 | shuttered | 12/12 shuttered | 12/12 shuttered | trail 11s -0.09; sky 117s 0.67 | shuttered | flats | unambiguous | – |
| Europa | 20251014 | shuttered | 7/7 shuttered | 7/7 shuttered | sky 61s 0.88; sky 108s 0.53 | shuttered | flats | unambiguous | – |
| Europa | 20251030 | shuttered | 7/7 shuttered | 7/7 shuttered | trail 13s -0.49 sh; trail 27s 2.79 | shuttered | flats | unambiguous | – |
| Europa | 20260106 | shuttered | 7/7 shuttered | 7/7 shuttered | trail 27s 0.34; sky 66s -0.20 sh | shuttered | flats | unambiguous | – |
| Europa | 20260110 | shuttered | 12/12 shuttered | 12/12 shuttered | sky 48s -0.82 sh; sky 142s 0.01 | shuttered | flats | unambiguous | – |
| Europa | 20260208 | shuttered | – | 11/16 shuttered | sky 66s 0.61; sky 80s -0.31 | shuttered | flats | probable | – |
| Europa | 20260316 | shuttered | 20/23 shuttered | 14/18 shuttered | trail 39s 0.25; sky 143s 0.20 | shuttered | flats | unambiguous | – |
| Europa | 20260505 | shuttered | 14/22 shuttered | – | trail 10s 0.04; trail 34s 0.35; sky 56s -0.41 | shuttered | flats | probable | – |
| Europa | 20260907 | shuttered | 0/1 inc. | – | sky 100s -0.19 | unknown | – | inconclusive | – |
| Ganymede | 20190118 | shuttered | 7/7 shuttered | – | trail 15s -0.85 sh; trail 25s 7.23 | shuttered | flats | probable | – |
| Ganymede | 20210401 | shuttered | 7/7 shuttered | – | sky 69s 0.60 | shuttered | flats | probable | – |
| Ganymede | 20230130 | shuttered | 7/7 shuttered | 7/7 shuttered | trail 10s 0.20; trail 31s -5.11 | shuttered | flats | unambiguous | – |
| Ganymede | 20241111 | shuttered | – | 7/7 shuttered | trail 30s 1.13; sky 51s 0.61 | shuttered | flats | probable | – |
| Ganymede | 20250323 | shuttered | – | – | trail 24s 0.92 | unknown | – | inconclusive | – |
| Ganymede | 20250404 | shuttered | 0/1 inc. | 13/18 shuttered | trail 10s 0.45; trail 20s 0.09 | shuttered | flats+science | unambiguous | – |
| Ganymede | 20250405 | shuttered | 0/1 inc. | 13/18 shuttered | trail 10s 0.35; trail 20s 0.23 | shuttered | flats+science | unambiguous | – |
| Ganymede | 20250406 | shuttered | – | 13/18 shuttered | trail 10s 0.55; trail 20s 0.27 | shuttered | flats+science | unambiguous | – |
| Ganymede | 20250407 | shuttered | – | 13/18 shuttered | trail 10s 0.61; trail 20s 0.25 | shuttered | flats+science | unambiguous | – |
| Ganymede | 20250408 | shutterless 1MHz | – | – | trail 20s 4.58 SL | shutterless 1MHz | science | probable | – |
| Ganymede | 20250409 | shutterless 1MHz | – | – | trail 10s 5.15 SL; trail 20s 4.32 SL | shutterless 1MHz | science | probable | – |
| Ganymede | 20250410 | shuttered | – | 11/16 shuttered | trail 20s 0.11; sky 72s -3.64 | shuttered | flats+science | unambiguous | – |
| Ganymede | 20250411 | shuttered | – | – | trail 20s 0.17 | unknown | – | inconclusive | – |
| Ganymede | 20250412 | shuttered | – | 14/18 shuttered | trail 20s 0.14; sky 42s -0.07 sh | shuttered | flats+science | unambiguous | – |
| Ganymede | 20250510 | shuttered | 20/20 shuttered | 16/18 shuttered | sky 41s 2.05 SL; sky 62s 0.79 | shuttered | flats | unambiguous | – |
| Ganymede | 20250608 | shuttered | 20/20 shuttered | 12/20 shuttered | sky 47s 0.62; sky 99s 1.09 | shuttered | flats | unambiguous | – |
| Ganymede | 20250829 | shuttered | – | – | – | unknown | – | inconclusive | – |
| Ganymede | 20251005 | shuttered | 19/20 shuttered | 10/15 shuttered | trail 36s 1.09; sky 117s 2.67 | shuttered | flats | unambiguous | – |
| Ganymede | 20251128 | shuttered | – | – | sky 46s 0.74; sky 156s 2.12 | unknown | – | inconclusive | – |
| Ganymede | 20251221 | shuttered | 24/30 shuttered | – | – | shuttered | flats | probable | – |
| Ganymede | 20260210 | shuttered | 21/24 shuttered | 13/18 shuttered | sky 66s 1.18; sky 140s 2.14 | shuttered | flats | unambiguous | – |
| Ganymede | 20260211 | shuttered | – | 9/13 shuttered | sky 66s 1.25; sky 140s 2.02 | shuttered | flats | probable | – |
| Ganymede | 20260212 | shuttered | – | 10/14 shuttered | sky 66s 1.25; sky 140s 2.10 | shuttered | flats | probable | – |
| Ganymede | 20260213 | shuttered | – | 13/17 shuttered | sky 140s 2.03 | shuttered | flats | probable | – |
| Ganymede | 20260214 | shuttered | – | 11/16 shuttered | sky 66s 1.08; sky 140s 1.93 | shuttered | flats | probable | – |
| Ganymede | 20260215 | shutterless 1MHz | – | – | – | unknown | – | inconclusive | – |
| Ganymede | 20260216 | shutterless 1MHz | – | – | trail 39s 5.42 SL | shutterless 1MHz | science | probable | – |
| Ganymede | 20260217 | shutterless 1MHz | – | – | – | unknown | – | inconclusive | – |
| Ganymede | 20260218 | shutterless 1MHz | 2/2 shutterless 1MHz (4.42) | 13/14 shutterless 1MHz (4.35) | trail 39s 5.19 SL | shutterless 1MHz | flats | unambiguous | – |
| Ganymede | 20260219 | shutterless 1MHz | – | 9/11 shutterless 1MHz (4.44) | trail 39s 5.13 SL | shutterless 1MHz | flats+science | unambiguous | – |
| Ganymede | 20260220 | shutterless 1MHz | 22/22 shutterless 1MHz (4.28) | 12/14 shutterless 1MHz (4.45) | trail 39s 5.05 SL | shutterless 1MHz | flats | unambiguous | – |
| Ganymede | 20260221 | shutterless 1MHz | 21/21 shutterless 1MHz (4.28) | 11/13 shutterless 1MHz (4.44) | trail 39s 5.34 SL | shutterless 1MHz | flats | unambiguous | – |
| Ganymede | 20260222 | shutterless 1MHz | – | – | sky 41s 5.40 SL | shutterless 1MHz | sky_ramp | probable | – |
| Ganymede | 20260223 | shutterless 1MHz | – | – | sky 41s 5.33 SL | shutterless 1MHz | sky_ramp | probable | – |
| Ganymede | 20260224 | shutterless 1MHz | 23/23 shutterless 1MHz (4.28) | 14/16 shutterless 1MHz (4.45) | sky 41s 5.19 SL | shutterless 1MHz | flats | unambiguous | – |
| Ganymede | 20260225 | shutterless 1MHz | 24/24 shutterless 1MHz (4.28) | 14/15 shutterless 1MHz (4.46) | sky 41s – | shutterless 1MHz | flats | unambiguous | – |
| Ganymede | 20260226 | shutterless 1MHz | 23/23 shutterless 1MHz (4.28) | 13/14 shutterless 1MHz (4.41) | sky 41s – | shutterless 1MHz | flats | unambiguous | – |
| Ganymede | 20260227 | shutterless 3MHz | 17/27 shutterless 3MHz (1.47) | 5/5 shutterless 3MHz (1.59) | sky 41s – | shutterless 3MHz | flats | unambiguous | – |
| Ganymede | 20260228 | shutterless 3MHz | 18/28 shutterless 3MHz (1.47) | – | sky 41s – | shutterless 3MHz | flats | probable | – |
| Ganymede | 20260301 | shutterless 3MHz | 18/28 shutterless 3MHz (1.47) | 12/19 shutterless 3MHz (1.50) | sky 140s –; sky 161s – | shutterless 3MHz | flats | unambiguous | – |
| Ganymede | 20260302 | shuttered | 20/30 shuttered | 0/1 inc. | trail 37s 1.35; sky 43s – | shuttered | flats | probable | – |
| Ganymede | 20260303 | shuttered | 20/30 shuttered | 0/1 inc. | trail 37s 1.38; sky 43s 0.02 | shuttered | flats | probable | – |
| Ganymede | 20260304 | shuttered | 20/23 shuttered | 13/18 shuttered | trail 25s -0.03; sky 157s 2.64 | shuttered | flats | unambiguous | – |
| Ganymede | 20260305 | shuttered | 20/23 shuttered | 10/15 shuttered | trail 37s 1.17; sky 157s -0.94 | shuttered | flats | unambiguous | – |
| Ganymede | 20260306 | shuttered | 21/24 shuttered | – | trail 25s 0.09 | shuttered | flats+science | unambiguous | – |
| Ganymede | 20260308 | shuttered | 20/23 shuttered | 21/24 mixed | trail 37s 0.72; sky 140s 5.26 | shuttered | flats+science | unambiguous | – |
| Ganymede | 20260611 | shuttered | – | 14/19 shuttered | sky 100s 0.68 | shuttered | flats | probable | – |
| Ganymede | 20260701 | shuttered | 15/15 shuttered | 10/17 shuttered | trail 20s -0.32; sky 66s 1.05 SL; sky 134s 1.60 | shuttered | flats | unambiguous | – |
| Ganymede | 20260707 | shuttered | 15/15 shuttered | 10/17 shuttered | trail 20s 2.14; sky 62s 1.09; sky 66s 1.76 SL | shuttered | flats | unambiguous | – |
| Ganymede | 20260827 | shuttered | – | – | sky 74s 0.41; sky 97s 1.25 | unknown | – | inconclusive | – |
| Ganymede | 20261001 | shuttered | 14/18 shuttered | 7/12 shuttered | sky 115s 1.78; sky 156s 1.88 | shuttered | flats | unambiguous | – |
| Io | 20180716 | shuttered | – | 15/15 shuttered | trail 23s 0.64; sky 50s 0.20 sh | shuttered | flats | probable | – |
| Io | 20180724 | shuttered | 10/10 shuttered | 10/10 shuttered | trail 20s -0.94 sh; trail 23s 0.31 | shuttered | flats | unambiguous | – |
| Io | 20210919 | shuttered | 10/10 shuttered | 10/10 shuttered | sky 96s 2.74; sky 180s – | shuttered | flats | unambiguous | – |
| Io | 20220406 | shuttered | – | 10/10 shuttered | trail 20s -0.98 sh | shuttered | flats | probable | – |
| Io | 20230222 | shuttered | 12/12 shuttered | – | – | shuttered | flats | probable | – |
| Io | 20250107 | shuttered | – | 7/7 shuttered | trail 10s -0.15; sky 61s -1.08 sh | shuttered | flats+science | unambiguous | – |
| Io | 20250202 | shuttered | 12/12 shuttered | 12/12 shuttered | trail 10s 0.22; trail 31s -0.28 | shuttered | flats | unambiguous | – |
| Io | 20250520 | shuttered | 12/12 shuttered | 12/12 shuttered | trail 23s -1.33 sh; sky 71s 0.80 | shuttered | flats | unambiguous | – |
| Io | 20250616 | shuttered | 17/17 shuttered | – | trail 10s 0.57; sky 90s 0.91 | shuttered | flats | probable | – |
| Io | 20250903 | shuttered | – | – | trail 30s 0.96 | unknown | – | inconclusive | – |
| Io | 20251029 | shuttered | 12/12 shuttered | 12/12 shuttered | trail 16s 0.22; trail 20s 1.22; trail 40s 0.17 | shuttered | flats | unambiguous | – |
| Io | 20251124 | shuttered | 17/17 shuttered | 17/17 shuttered | trail 14s -0.14; trail 24s -0.15; sky 90s 0.83 | shuttered | flats | unambiguous | – |
| Io | 20251126 | shuttered | 17/17 shuttered | 17/17 shuttered | trail 14s 0.17; trail 24s -0.20 sh; sky 100s 1.73 | shuttered | flats | unambiguous | – |
| Io | 20251225 | shuttered | 6/7 shuttered | 5/7 shuttered | sky 45s 1.03 SL; sky 170s 0.25 | shuttered | flats | unambiguous | – |
| Io | 20260121 | shuttered | 7/7 shuttered | 4/7 shuttered | trail 28s -1.97; sky 177s 2.92 | shuttered | flats | unambiguous | – |
| Io | 20260424 | shuttered | 22/29 shuttered | 12/17 shuttered | trail 23s 2.83; sky 159s 22.39 | shuttered | flats | unambiguous | – |
| Io | 20260525 | shuttered | 20/24 shuttered | – | sky 129s 5.68; sky 157s 2.23 | shuttered | flats | probable | – |


