# The WCS cache: T7's plate solutions, kept outside the raw tree

T7 plate-solves every calibrated frame (`astrom/astrometry.py`, `astrom/pointer_wcs.py`). It keeps each night's
solutions in a table of its own, so that a rerun of the night can apply them instead of solving again. T7 never
writes a raw frame.

## Why

Before the cache, T7 copied every solution back into the raw frame (`copy_wcs_to_raw`), so that the next run would
find a WCS in the frame's header and skip it. That cache did not work:

- **SSO never used it.** T7's "already solved" test could not tell Astra's on-site WCS from orchard's, and Astra
  frames arrive with one, so SSO always ran with `--force-platesolve`. In 30 days of logs, about 379,000 raw-frame
  writes were never read. Only Artemis skipped frames, about 3,800 a month.
- **It edited raw frames with no record.** The copy merged orchard's cards into Astra's and kept Astra's
  comments.
- **Any tree could write it.** A scratch tree pointed at the production data root wrote into production raw frames.

## Where

One table per night:

```
<BASEDIR>/PipelineOutput/wcs_cache/<TEL>/<YYYYMMDD>.ecsv
```

- ZLP works the path out in T7 and passes it as `astrometry.py --wcs-cache`. Set `WCS_CACHE_DIR` to use
  `$WCS_CACHE_DIR/<TEL>/<YYYYMMDD>.ecsv` instead. ZLP's own command line is unchanged.
- The table sits outside `v2`/`v3`, so reruns, T12 and promotion never move it. A scratch data root has its own.
- `astrometry.py` without `--wcs-cache` solves every frame and keeps nothing.

## What a row holds

One row per raw frame, keyed by the raw file's name. Calibrated frames are matched to raw frames by T6's naming
rule (`pipeutils.proc_name`), through T1's image list (`--raw-list`, read only).

| Columns | |
|---|---|
| `raw` | the raw file's name |
| `raw_size`, `raw_mtime_ns` | the raw file when it was solved. A re-download, rotation or repair changes them. |
| `naxis1`, `naxis2`, `trim_x`, `trim_y` | the calibrated frame's shape, and where it sits in the raw frame |
| `solver_epoch` | `astrometry.SOLVER_EPOCH` (see below) |
| `gaia_db`, `gaia_db_size`, `gaia_db_mtime_ns` | the Gaia database it was solved against |
| `astropy`, `twirl` | the library versions |
| `orchver`, `orchrel`, `commit` | the orchard that solved it: version, code release (empty outside a release), commit |
| `solved_at`, `matches`, `proc` | when (UTC), how many stars matched, which calibrated frame |
| `cards` | the solution: every card the solve wrote (the WCS keywords and `MJDREF`), 80-character card images exactly as they stand in the solved frame |

Only successful solves get a row. A failed or timed-out frame is solved again on the next run.

## When a row is used

All of these must hold:

1. Every validity field (`raw_size` to `twirl` above) matches the frame and the run.
2. A code release wrote the row (`orchrel` is set), or the commit that is running now wrote it. A scratch tree
   (staged with `git archive`, so no release) therefore never feeds a production run, even on the production data
   root. Its own rows feed only itself. Rows from a checkout with uncommitted changes (`-dirty`) are never used.
3. The cards parse as a complete celestial WCS.

`--force-platesolve` ignores the table: every frame is solved and every solved frame's row is replaced.

A WCS already in a frame's header (Astra's, PinPoint's, ESO's, or orchard's own in raw frames written by older
code) never counts as a solution. Without a usable row, a frame is always solved.

## What T7 does

1. Reads the night's table once.
2. For each frame with a usable row, applies it, header only: it clears the WCS keywords as a solve does, applies
   the cards with `header.update()`, and records the edit (`ORCHMOD`, HISTORY `orchard <ver> T7: WCS from cache`). The
   pixels are never read, so a 16-bit frame stays 16-bit. The result is the header a fresh solve gives; only that
   HISTORY card differs from a solve's `T7: WCS solved`.
3. Solves the other frames in the pool, as before.
4. After the pool, the parent reads the solved frames' headers and writes their rows: one read-change-write under
   `utils/shared_files.locked`, through a temporary file and `os.replace`. Two runs on one night (two targets, a
   rerun beside the nightly run) lose no rows, and a reader sees the old table or the new one.

Nothing here can end a night: T7 runs under `errexit` in ZLP. An unreadable table, a bad row, a missing Gaia
database or an unwritable directory is a `WARNING: WCS cache ...` line, and the frames are solved. An unreadable
table is moved to `<table>.corrupt` and replaced by the next write.

The log says what happened:

```
WCS cache: .../wcs_cache/Io/20261018.ecsv (1012 rows)
WCS cache: 1009 frames from the cache, 3 to solve (no entry 3)
...
WCS cache: 2 rows written to .../wcs_cache/Io/20261018.ecsv
```

The reasons in brackets are `no entry`, the first validity field that differs (e.g. `raw_mtime_ns`),
`provenance`, `forced`, `no raw frame` (not in T1's list) or `unusable` (the row did not apply).

## SOLVER_EPOCH

`SOLVER_EPOCH` in `astrom/astrometry.py` is an integer recorded in every row. Bump it in any PR that changes the
solutions T7 writes: `pointer_wcs`'s detection, catalogue query, matching or acceptance, or the twirl call. Rows
from the old solver then stop being used, and each night is solved afresh on its next run. It follows the pattern
of `readout.RECORD_VERSION`.

A change to the Gaia database, astropy or twirl invalidates rows by itself.

## Looking at a table

```python
from astropy.table import Table
t = Table.read('/data/SPECULOOSPipeline/PipelineOutput/wcs_cache/Io/20261018.ecsv', format='ascii.ecsv')
t['raw', 'orchver', 'orchrel', 'solved_at', 'matches']
```

`astrom/wcs_cache.py` has `load()`, which returns `{raw name: row}`, and `cards_header()`, which turns a row's
cards into an astropy `Header`.

Deleting a night's table is always safe: the next run solves the night and writes it again.

## Rollout

1. The cache deploys with SSO still forced (`force_platesolve` is true for Io, Europa, Ganymede and Callisto in
   `jobqueue/config.py`; the production `queue/config.json` does not override it). The only changes are that raw
   frames are no longer written and every run fills the table.
2. After validation, a small follow-up sets `force_platesolve` to false for the SSO telescopes, so reruns take
   their solutions from the table. Solutions can differ slightly from a re-solve, because the Gaia shard order in
   `pointer_wcs.db_query` is not deterministic.

The first run of each Artemis night after that deploy solves every frame, including frames whose raw files carry
orchard's WCS from the old copy. That happens once per night.

Raw frames already written by the old code keep the cards it wrote; nothing removes them.
