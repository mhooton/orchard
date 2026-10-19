# Versioning: which orchard made a product

Orchard's versions are git tags on `main`, `vMAJOR.MINOR.PATCH`. Every data product the pipeline builds a header
for carries the version of the code that wrote it (`ORCHVER`), and every queue job also carries the code release it
ran from (`ORCHREL`). So for any product you can answer "does this predate fix X?" or "which products must be
reprocessed?" from the file alone.

## The version string

A piece of code's version is `git describe` of the commit it came from, without the leading `v`:

| `ORCHVER` | Meaning |
|---|---|
| `3.1.0` | The tagged commit itself. |
| `3.1.0-4-g9a0de6d` | Four commits after `v3.1.0`, at commit `9a0de6d`. Exact, but not a numbered release. |
| `9a0de6d` | No version tag is reachable from the commit. |
| `...-dirty` | Run from a checkout with uncommitted changes. Never seen in production. |
| `unknown` | The code was neither a git archive nor a checkout (see "How the version reaches production" below). |

A product made before versioning has no `ORCHVER` at all.

`VERSION`, `PV` and `PIPE_V` (`'v3'`, or a hardcoded `'v2'` in the PWV files) are a different thing: they name the output tree
(`PipelineOutput/v3`), not the code. They stay as they are, because the promotion checker reads `PIPE_V`.

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
git tag -a v3.2.0 -m "orchard 3.2.0: stable ALC weights; stale PC matrix fix" <sha>
git push origin v3.2.0
```

A tag made on GitHub's release page counts too; lightweight tags are fine. A commit deployed without a tag is
stamped `3.1.0-N-g<sha>`. That still identifies it exactly; it just has no number of its own.

## Where it appears

**FITS keywords**, in the primary header:

| Product | Stage | Keywords |
|---|---|---|
| `<run>_MasterBias.fits` | T2 | `ORCHVER`, `ORCHREL` |
| `<run>_MasterDark*.fits` | T3 | `ORCHVER`, `ORCHREL` |
| `<run>_MasterFlat_<filt>.fits` | T4 | `ORCHVER`, `ORCHREL` |
| `<run>_BadPixelMap.fits`, `hot/cold/flatbad_pixel_map.fits` | T5 | `ORCHVER`, `ORCHREL` |
| `proc*.fits` (calibrated frames) | T6 | `ORCHVER`, `ORCHREL` (later stages edit these frames in place and keep them) |
| `<dr2>_<filt>_<night>_output.fits` | T9 | `ORCHVER`, `ORCHREL` |
| `<dr2>_<filt>_output.fits` (global, every night) | T9 | Version of the **last** run that updated it, with that comment, as `VERSION` already behaves |
| `<gaia>_<filt>_<night>_<ap>_diff.fits` | T10 | `ORCHVER`, `ORCHREL` |
| `<gaia>_<filt>_<night>_pwv_*.fits` | T10 | `ORCHVER`, `ORCHREL` |
| `StackImages_nightly/<night>/<target>_stack.fits` | T11 | `ORCHVER`, `ORCHREL` |

`ORCHREL` is written only by code running from a release. Code run from the live tree outside the queue, or from a
scratch tree, writes `ORCHVER` alone.

To read the keywords:

```bash
python -c "from astropy.io import fits; import sys; print(fits.getval(sys.argv[1], 'ORCHVER'))" <file>
```

**Elsewhere:**

- The PDF night report (T11): the version line is under the title on page 1, and in the PDF's Creator field.
- The pipeline log: `orchard 3.1.0 (commit 9a0de6d..., release 20261019T...)` near the top of every run.
- The job queue: `q status` and `q release list` show each release's version, and so does `RELEASE.json`.

**Not stamped**, deliberately:

- Files CASU writes: stack catalogues, `.cat`, `.phot`, stack images. Their headers come from the tool or its
  inputs.
- Text products. `*_bestap.txt`, `readoutnoise.dat`, `darkcurrent.dat` and the `_MCMC` light curves are parsed
  by position, by the pipeline and by readers outside it, and a comment line would break them.

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
python src/utils/orchard_version.py            # 3.1.0-4-g9a0de6d
python src/utils/orchard_version.py --line     # orchard 3.1.0-4-g9a0de6d (commit 9a0de6d..., release ...)
```
