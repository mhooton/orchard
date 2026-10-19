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
| `proc*.fits` (calibrated frames) | T6 | yes | T7 `inherited WCS removed`, `WCS solved`; T8/T9 `wcsfit: WCS refit` (or `... failed`); framepack (`ORCHMOD` only; its own `framepack:` card is the record) |
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
