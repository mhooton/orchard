# Plate solving (T7)

T7 gives every calibrated science frame a WCS (`astrom/astrometry.py`, which runs `astrom/pointer_wcs.py` on each
frame). A frame without one never reaches the light curve: T9's CASU `wcsfit` starts from it, and T10 drops the
frame. So a frame T7 cannot solve is a missing point, and a frame it solves wrongly is a point measured on the wrong
pixels.

## How a frame is solved

1. **Detections.** `clean_image` and `find_stars_multiscale`, as before: brightest first.
2. **Catalogue.** One query of the local Gaia/2MASS database, in a box twice the field about the header pointing
   (`known_scale.CATALOGUE_BOX`). The stars inside the old 1.2-field box are what the log's `Gaia: N` still counts,
   and what twirl gets if it runs.
3. **The frame's inherited WCS**, if it has one (Astra's own solve, or an earlier T7 solution copied into the raw
   frame), is recentred by an offset vote, refined, and kept only if it passes the verification below. With
   `--force-platesolve` it is read before T7 clears it.
4. **Known-scale solve** (`astrom/known_scale.py`). The plate scale is known to about 1% (header, or instrument config)
   and the pixels are square, so the sky-to-pixel map is a rotation, a parity and a shift:
   - **pairs**: a detection pair and a catalogue pair of the same length vote for the angle between them. Each
     candidate rotation (both parities, both directions) is scored by voting on the offset between every detection
     and every catalogue star; the true one gives a sharp peak.
   - **grid**: in crowded fields the pair votes drown, so a 1-degree rotation grid, scored the same way, is the
     fallback (about 5 s on a bulge field).
   - The winner is refined by iterative cross-matching: a TAN fit with 10 or more matched stars, a similarity
     (scale, rotation, shift) below that.
5. **twirl** (the solver until 2026-10), unchanged, only when the known-scale solve has nothing verified.
6. **A second pass** (`astrometry.py`), after every frame has had its turn. A frame that failed with 5 or more round
   sources, not trailed and not timed out, is tried again from the solution of the frame nearest to it in time that
   has one: solved this run, or applied from the WCS cache, and no more than 2 hours away. The neighbour's WCS is
   recentred within 50 px and refined (`known_scale.solve_from_neighbour`; with too few matches to refit, it is kept
   as recentred). The result is kept if it verifies as below, or if the sources are round and 3 or more sit on
   catalogue stars at chance 10^-8 or less and the solution agrees with the neighbour's (centre within 50 px,
   rotation within 0.5 degrees, scale within 1%). Faint frames of sparse fields have only 5-7 sources, and alone they
   cannot verify. On Callisto 20261012 (SPIRIT) v4.1.0 refused 63 such frames whose old photometry was normal; their
   consecutive frames, 6-18 s apart, put 3 of their 5 sources on catalogue stars (twirl's "4 of 5" was its own fitted
   quad). A frame the second pass cannot solve goes through steps 3-5 again and keeps its first failure message.
   Solved lines say `Method: neighbour`, and the solutions go into the WCS cache like any other (`docs/wcs-cache.md`).

## The verification

Every solution, whichever step made it, must pass `known_scale.verified()`. Of the brightest 100 detections:

- at least 6 lie within 3 px of a catalogue star,
- that many is at most 10^-10 likely by chance, for the catalogue density in the frame,
- and they are at least a quarter of what could match (`min(detections tested, catalogue stars in the frame)`), or
  15% when the sources are round (`trail_ellipticity` below 0.5).

In a sparse field 5 matches also verify, if they are at least 80% of what could match. Four never do.

A known-scale solution must also be what a telescope gives: the expected scale to 5% and, for a TAN fit, square
pixels to 2%.

Why these numbers:

- **Six, or five nearly complete.** Two matched stars fix a similarity, so only the matches beyond two are
  evidence. In a frame holding 7 catalogue stars, two chance coincidences beyond the first two turn up about once in
  ten frames, which is why 4 matches is never enough. Three more coincidences, among only 5 or 6 stars, are
  negligible. Old 5-match solutions in sparse SPIRIT fields agreed with this solve's.
- **The share of matches** is there for star trails. A trail holds each star's position at some instant, so on frames
  from a telescope that stopped tracking the best wrong candidates reached 8-19% of what could match, at chance
  10^-10 to 10^-21. Correct solutions match 70-100% on most nights, and down to 21% on Artemis I+z frames, where most
  bright detections are fainter than the catalogue's J = 16. Those frames' sources are round, so 15% applies.

## Reading a failure

The per-frame line in the log keeps its form; a solved line now names the method:

```
[12/915] procX.fits SUCCESS - Sources: 218→100 Gaia: 75→61 Matches: 48 Method: pairs CRPIX=(...) CRVAL=(...)
[13/915] procY.fits FAILED - Sources: 372→16 Gaia: 270→32 Matches: 0 Error: <why>
```

For a known-scale solution `Sources: a→b` is detections → detections tested, `Gaia: c→d` is the 1.2-field box →
catalogue stars in the frame, and `Matches` is how many tested detections sit on catalogue stars.

| Error starts with | Meaning | Whose loss |
|---|---|---|
| `Not enough stars detected` | 4 or fewer sources: cloud, or the dome | the sky |
| `Sources trailed` | the brightest sources all stretched one way: the telescope stopped tracking | the sky |
| `No verified solution (best m of n ...)` | nothing verified; twirl's own reason follows | the pipeline, unless the field is too sparse (see below) |
| `Not enough Gaia stars found in field` | 4 or fewer catalogue stars in the box | the catalogue |

`utils/promote_policy.py` counts the first two as the sky's losses (`DATA_ERRORS`) and the rest, on frames with 20
or more sources, as the pipeline's.

## What it does not fix

- **Fields too sparse for the catalogue.** The local database stops at J = 16. A SPIRIT field with 4 catalogue
  stars in the frame (Callisto 20250810 Sp0055-3052) cannot reach 6 matches, and 4 is not evidence.
- **The telescope pointing elsewhere.** On Artemis 20240211 the field walked 4, 6.5 and then 15 arcmin from the
  target over the night; frames more than a field-width from the header pointing are outside the catalogue box, and
  the target is outside them.
- **Star trails.** Frames from a telescope that stopped tracking are reported, not solved.

## One BLAS thread per process

T7 runs its frames in a process pool (20 processes in production), and `astrometry.py` sets `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS` and `MKL_NUM_THREADS` to 1 before numpy is first imported. Keep that assignment above every other
import. Without it, a forked worker whose OpenBLAS starts its thread pool can deadlock. From v4.1.0 to v4.3.1,
10 of 29 T7 runs on appct stalled that way (2026-10-20/21): workers with 64 threads, all asleep in `futex_wait`, which
the 60 s per-frame alarm cannot interrupt, so the job sat until the queue's timeout. numpy and scipy each load their
own OpenBLAS, and both read these variables.

## Measured

See the PR that introduced this (2026-10-19) for the before/after on 24 nights; the scratch work is in
`mh_scratch/platesolve_partial/` on appct.
