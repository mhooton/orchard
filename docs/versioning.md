# Versioning: which orchard made a product

Orchard's versions are git tags on `main`, `vMAJOR.MINOR.PATCH`, starting at `v4.0.0`. That number separates the
versioned code from the loosely defined "v3" before it. The `v2` and `v3` folders under `PipelineOutput` are only
the mechanism that controls what the portal ingests; they say nothing about the code.

Every FITS file a pipeline stage creates records the version that made it (`ORCHVER`). Every later in-place edit
of that file records itself too (`ORCHMOD` and a HISTORY card). So from the file alone you can answer "does this
predate fix X?" and "which products must be reprocessed?".

## The version string

A piece of code's version is `git describe` of the commit it came from, without the leading `v`:

| Version | Meaning |
|---|---|
| `4.1.0` | The tagged commit itself. |
| `4.1.0-4-g9a0de6d` | Four commits after `v4.1.0`, at commit `9a0de6d`. Exact, but not a numbered release. |
| `9a0de6d` | No version tag is reachable from the commit. |
| `...-dirty` | Run from a checkout with uncommitted changes. Never seen in production. |
| `unknown` | The code was neither a git archive nor a checkout (see "How the version reaches production" below). |

A file made before versioning has no `ORCHVER` at all.

`VERSION`, `PV` and `PIPE_V` (`'v3'`, or a hardcoded `'v2'` in the PWV files) are a different thing: they name the
output tree, not the code. They stay as they are, because the promotion checker reads `PIPE_V`.

## When the number goes up

The question each part answers is "do my products change?":

- **MAJOR**: a product's format or meaning changes in a way that breaks readers. Examples are renamed or removed
  columns or extensions, changed units, and a changed directory layout.
- **MINOR**: the values in products can change for nights that already processed. Examples are a new correction,
  a calibration fix, and changes to the ALC, the apertures or the target identification. Products from before and
  after a minor bump are not comparable without checking.
- **PATCH**: products are unchanged for nights that already worked. Examples are crash fixes, queue and operations
  changes, reporting, and docs.

PRs do not touch the version. A PR that bumps it would conflict with every other open PR. The bump happens when
something is deployed: tag the commit being deployed, then stage it. Several merged PRs can share one tag.

```bash
git tag -a v4.1.0 -m "orchard 4.1.0: <what changed in the products>" <sha>
git push origin v4.1.0
```

A tag made on GitHub's release page counts too; lightweight tags are fine. A commit deployed without a tag is
stamped `4.0.0-N-g<sha>`. That still identifies it exactly; it just has no number of its own.

## The keywords

All of them are in the primary header (HDU 0), or HDU 1 of a Rice-compressed `.fits.fz` frame.

| Keyword | Set | Meaning |
|---|---|---|
| `ORCHVER` | When the file is created; never changed afterwards | The orchard version that made the file. |
| `ORCHREL` | When the file is created, by code running from a queue release | The code release that made it, e.g. `20261019T112705Z-d75c5a1`. |
| `ORCHMOD` | By every later in-place edit | The orchard version of the latest edit. |
| `HISTORY orchard <version> <script>: <what>` | By every later in-place edit | One card per edit. An identical card is not added twice, so T8 and T9 refitting a frame's WCS leave one. |

A file a tool makes by copying another file's header is stamped as new. This covers CASU's stack catalogue,
`outstack`, `.cat` and `.phot`, and the T11 stack, which copies a proc frame's header. Such a file gets its own
`ORCHVER`, and the copied `ORCHMOD` and `ORCHREL` are dropped. The copied HISTORY cards stay, as the history of
the frame it was made from.

The edit cards never repeat text that code searches HISTORY for: `horizontal mirroring`,
`anticlockwise rotation`, `orient repair` and `framepack:`. Each card fits in one card (72 characters). The
helpers live in `src/utils/orchard_version.py`:

- `stamp` and `stamp_file` for creation.
- `record_edit`, `edit_file` and `edit_fitsio_hdu` for edits.

The file-level ones never raise: a stamp that fails is logged, and the stage carries on.

## Where it appears

| File | Stage | Created (`ORCHVER`) | Later edits (`ORCHMOD` + HISTORY) |
|---|---|---|---|
| `<run>_MasterBias.fits`, `_MasterDark*.fits`, `_MasterFlat_<filt>.fits` | T2–T4 | yes | — |
| `<run>_BadPixelMap.fits`, `hot/cold/flatbad_pixel_map.fits` | T5 | yes | — |
| `proc*.fits` (calibrated frames) | T6 | yes | T7 `inherited WCS removed`, `WCS solved` or `WCS from cache` ([wcs-cache.md](wcs-cache.md)); T8/T9 `wcsfit: WCS refit` (or `... failed`); framepack (`ORCHMOD` only; its own `framepack:` card is the record) |
| `<target>/1_initial-catalogue.fits`, `proc*.cat` | T8/T9 | yes | — |
| `*_outstack_<filt>.fits` (+ `_conf`), `*_stack_catalogue_<filt>.fits` | T8 | yes | A backup catalogue restored for a later night: `catalogue_fov: target roles for <night>`. Fallback crossmatch: `proper_motion: Gaia crossmatch` |
| `*_stack_catalogue_pm_<filt>.fits` | T9 | yes | — |
| `proc*.fits.phot` | T9 | yes | — |
| `<dr2>_<filt>_<night>_output.fits` | T9 | yes | — |
| `<dr2>_<filt>_output.fits` (global, every night) | T9 | by the run that made it | `ORCHMOD` only, on every update: HISTORY there is the list of nights, and `SPlightcurve` reads `HISTORY[0]` as a date |
| `<gaia>_<filt>_<night>_<ap>_diff.fits`, `..._pwv_*.fits` | T10 | yes | — |
| `StackImages_nightly/<night>/<target>_stack.fits` | T11 | yes | — |

To read the keywords:

```bash
python -c "from astropy.io import fits; import sys; h = fits.getheader(sys.argv[1]); print(h.get('ORCHVER'), h.get('ORCHMOD'))" <file>
```

**Elsewhere:**

- The PDF night report (T11): the version line is under the title on page 1, and in the PDF's Creator field.
- The pipeline log: `orchard 4.0.0 (commit 9a0de6d..., release 20261019T...)` near the top of every run.
- The job queue: `q status` and `q release list` show each release's version, and so does `RELEASE.json`.

**Not stamped**, deliberately: text products. `*_bestap.txt`, `readoutnoise.dat`, `darkcurrent.dat` and the
`_MCMC` light curves are parsed by position, by the pipeline and by readers outside it, and a comment line would
break them.

## Raw frames

Raw frames live in `Observations/<telescope>/images/<night>/`. ESO's downloads are kept exactly as delivered.
Every file orchard writes there records itself, by the same rules as products: `ORCHVER` when orchard writes the
frame at ingest, and `ORCHMOD` plus a HISTORY card for every later in-place edit.

| Frame | Written by | Keywords | HISTORY |
|---|---|---|---|
| Downloaded and needing no transform (ACP era, any non-Astra frame) | `request_eso.py`, the queue's fetch | none: byte for byte as ESO delivered it | — |
| Unpacked from a datacube (SPIRIT) | `unpack_datacubes.py`, at ingest | `ORCHVER` (+ `ORCHREL`) | `orchard <version> unpack_datacubes` |
| Astra frame rotated (Andor) or mirrored (SPIRIT) on arrival: lights and calibration frames | `astra_transform()`, called by `SSO_download.py` and the queue's fetch | `ORCHVER` (+ `ORCHREL`); `ASTRAROT`/`ASTRAMIR` and `ASTRAWCS` as before | the transform's own `WCS updated for ...` card (lights), then `orchard <version> astra_transform` |
| Transformed later from the command line | `python download/astra_transform.py ...` | `ORCHMOD` | `orchard <version> astra_transform: <mode> k=<n>` (or `: mirror`) |
| WCS repaired | `repair_astra_wcs.py --apply` | `ORCHMOD`; `ASTRAWCS = 2` as before | `repair_astra_wcs: <signatures>` as before, then `orchard <version> repair_astra_wcs: <signatures>` |
| Orientation repaired | `repair_orientation.py --apply` | `ORCHMOD` | `<date> orient repair: ...` as before, then `orchard <version> repair_orientation: rot90 k=<k>` (or `: fliplr`) |
| Packed or unpacked by hand | `utils/framepack.py` | `ORCHMOD` | its own `framepack: ...` card |

A frame that is both unpacked and turned on arrival carries `ORCHVER` once and both ingest cards, in that order.

How to read a raw frame's header:

- **`ORCHVER`**: orchard wrote this file at ingest, in that version.
- **`ORCHMOD`**: orchard has edited the file in place since. The `orchard ...` cards say which tools did it.
- **Neither**: the file is one of two things:
  - ESO's file as delivered;
  - a frame orchard wrote before raw frames were versioned. These carry `ASTRAROT`/`ASTRAMIR` and the
    descriptive HISTORY cards, and nothing more. The archive is not retrofitted.

One write is not recorded: T7's copy of its plate solution into raw light frames (`copy_wcs_to_raw`). It is
being replaced by a WCS cache outside the raw tree.

The descriptive cards stay exactly as they were, because the repair scripts read them back:

- `repair_orientation` counts the `horizontal mirroring` cards and refuses a frame that carries `orient repair`.
- `repair_astra_wcs` looks for `anticlockwise rotation` and `horizontal mirroring`.

The `orchard ...` cards never contain those phrases (`edit_card` raises if one would), so a stamped frame is
classified exactly as before.

### Checksums

Every orchard write to a raw frame leaves `CHECKSUM` and `DATASUM` valid:

- **Whole-file writers** pass `checksum=True` to astropy: `unpack_datacubes`, `astra_transform`,
  `repair_orientation` and framepack.
- **`repair_astra_wcs`** replaces header bytes only. It recomputes `CHECKSUM` from the new header and the frame's
  `DATASUM`:
  - If the old `CHECKSUM` agreed with that `DATASUM`, the `DATASUM` is trusted. This check reads the header only.
  - Otherwise `DATASUM` is computed from the data bytes, which are read but never decoded
    (`download/raw_provenance.py`).

Before this change, most of the archive's Astra frames did not verify. A survey on 2026-10-19 sampled 168 frames
over 15 nights on Io, Europa, Ganymede and Callisto
(`mh_scratch/raw-provenance/survey/` on appcs):

- **ESO's downloads verify**, ACP-era nights included.
- **Frames written by `astra_transform` don't**, on every night since Astra. `writeto(checksum=False)` writes back
  the cards it was given, so every transformed frame kept ESO's `CHECKSUM`/`DATASUM` for the image before it was
  turned. `repair_orientation` wrote the same way.
- **The exception is light frames that T7 later wrote a WCS into.** astropy's update mode recomputes the cards a
  file already has, so those frames verify, dated by T7's run.
- **Calibration frames on Astra nights never verify**, because T7 never writes them.
- **`repair_astra_wcs` header edits left `CHECKSUM` stale.**
- **Datacube frames verify only until the transform.** The unpacker's update-mode reopen recomputed the cube's
  copied cards, and the mirror then made them stale again.

So a raw frame's checksum can be relied on when the frame carries `ORCHVER` or `ORCHMOD`. Without them, it holds
only for ESO's untouched downloads. To check one (`1` valid, `0` wrong, `2` missing):

```bash
python -c "from astropy.io import fits; import sys; h = fits.open(sys.argv[1])[0]; print(h.verify_checksum(), h.verify_datasum())" <frame>
```

## How the version reaches production

Production runs from frozen copies of `src/`, not from a git checkout, so the code must carry its own version.
`src/utils/orchard_version.py` holds `$Format:...$` placeholders, and `.gitattributes` marks it `export-subst`.
`git archive` fills the placeholders with the commit's describe, hash and date, so every tree staged with
`git archive <sha> src` knows its version. That is the deploy route ("Deploying main" in
[job-queue.md](job-queue.md)), and scratch trees can be staged the same way. In a checkout the module asks `git describe` instead.

`q release build` and `q release deploy` read the staged tree's version and record it in `RELEASE.json`. They
refuse a stage whose archived commit is not the `--sha` given. A stage that was not made by `git archive` (for
example rsynced from a checkout) is built with a warning, and its products say `unknown` or the bare sha.

A git older than 2.35 cannot fill the describe placeholder. Products from such a stage get the bare sha, which is
still exact. Appcs and appct have git 2.6, so stage from a Mac or from the container (git 2.39).

```bash
python src/utils/orchard_version.py            # 4.0.0-4-g9a0de6d
python src/utils/orchard_version.py --line     # orchard 4.0.0-4-g9a0de6d (commit 9a0de6d..., release ...)
```
