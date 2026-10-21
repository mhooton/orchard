# Frame compression

Calibrated and raw frames are stored compressed, without losing anything the pipeline needs. Decisions D1-D9 of
the frame compression inventory (16 October 2026) and Matthew's choice of 21 October 2026 for raw frames.

| Frame | Stored as | Size | Exactness | Made by |
|---|---|---|---|---|
| Calibrated `<target>/1/procX.fits` | `procX.fits.fz`: 16-bit integers (per-frame BZERO, NaN as BLANK) + lossless RICE_1 | ~1/4 of float32 | pixels within 0.5 ADU of T6's float32 (light curves validated on four nights), then bit-exact | `pack_frames` after T11 (`--frame-format rice`, the default); the archive sweep |
| Raw `images/<night>/X.fits` / `X.fts` | `X.fits.fz` / `X.fts.fz`: lossless RICE_1 plus the original header block (`ORIGHDR`) | ~1/1.5-1/3 | the original file byte for byte: `ORIGSHA` (sha256), `ORIGSIZE`, `ORIGMTIM` checked on every decompression; ESO's `CHECKSUM` verifies on the result | `utils/framepack.py rawpack` (the raw sweep, not started) |

A compressed file is always the plain name plus `.fz` (D1), an empty primary HDU and the image in HDU 1 (as fpack
writes it). The 690 raw nights of 2019-20 that were compressed with fpack in January 2022 have no `ORIGHDR`: they
decompress to the same pixels and cards, not to the original bytes, which no longer exist.

## How the pipeline meets compressed frames

The pipeline works on plain frames (D2).

- **Raw frames.** When T1-T7 run, `ZLP_pipeline.sh` decompresses the night first (`utils/rawpack.py unpack`):
  - each working copy is the stored frame exactly, with its original mtime, so T7's WCS cache (keyed on the raw
    name, size and mtime) still matches;
  - the `.fz` waits in a hidden `.rawpack/` beside it, which no glob in the pipeline sees.
- **After the run.** When the date is done, and on any exit (error, the queue's SIGTERM, ^C), ZLP stores the night
  compressed again (`rawpack.py repack`):
  - a working copy that is unchanged (same inode, size, mtime and ctime, or the same bytes) is deleted and its
    `.fz` put back;
  - a working copy the run changed is compressed again, and the earlier `.fz` is kept as `.rawpack/X.fz.<stamp>`.
    Nothing in the pipeline writes raw frames since PR #64, so this should never happen.
- **A repack cut short** (the queue's SIGKILL two minutes after SIGTERM) leaves the night decompressed. Every
  reader handles that. The next run on the night takes over the working copies and finishes the repack.
- **Calibration frames from earlier nights.** T1 lists them whether or not they are compressed
  (`createlists.get_liste`), and T2-T6 read them through `pipeutils.open_fits_file`, which presents a `.fz` frame
  as the plain file in memory. Nothing outside the job's own night is decompressed on disk, so no other job's night
  is touched.
- **Calibrated frames.** T6 writes plain float32 frames, and removes an older run's `procX.fits.fz` beside one it
  writes. After T11, `pack_frames` converts them (`--frame-format rice`, the default; `int16` or `float32` to
  keep them plain).
- **Reruns.** A rerun from T7, T8 or T9 decompresses the target's packed frames first (`unpack_frames`), because
  those stages write WCS into them; `pack_frames` packs them again. A rerun from T10 or T11 reads them as they are:
  T10 follows its image list's plain names to the `.fz`, and T11 stacks `proc*.fits.fz`.
- **Restage and promotion.** `utils/stages.py` treats `procX.fits.fz` as the same product as `procX.fits`, so a
  packed night is restaged and judged like any other (S1, O3).

### Reading a frame that may be compressed

```python
from utils import framepack
with framepack.open_frame(path) as hdul:      # plain or .fz; a plain name that is now .fz is followed
    data, header = hdul[0].data, hdul[0].header
header = framepack.getheader(path)            # fits.getheader gives a .fz's empty primary header
```

`open_frame` never opens anything for writing. Never open a `.fz` with `mode='update'`: astropy rewrites a Rice
image it was asked to touch as lossy quantised float, without a warning.

## Commands

All of them run inside `orchard-server`, take the frames' directories, and verify every file before it replaces
the original (a hidden temporary file, read back and checked, then renamed; nothing on the servers is backed up).

```
python utils/framepack.py pack <target>/1/ --workers 8 --manifest <csv>      # calibrated: 16-bit + Rice
python utils/framepack.py rawpack <images/NIGHT> --workers 8 --manifest <csv> # raw: byte-exact Rice
python utils/framepack.py decompress <dir> --workers 8                       # either, back to plain
python utils/rawpack.py status|unpack|repack <images/NIGHT> --workers 8      # a raw night around a run
```

Converting a night takes it from the pipeline: hold the night (a queue job with the `night:` lock), never while a
pipeline run or promotion works on it.

## What is not compressed

- Raw frames that are not 16-bit integers: the 447 ESO float32 nights of Aug 2018-Feb 2019. `rawpack` refuses
  them, and they stay plain.
- Anything `rawpack` cannot rebuild byte for byte: extra HDUs or bytes after the image, non-zero padding.
- Master calibrations, stacks, catalogues, `StackImages_nightly` (D8).

## Known limits

- `download/repair_astra_wcs.py` rewrites header bytes in place, so it skips `.fz` frames and says so. Unpack a
  night (`rawpack.py unpack`) to repair it, then repack.
- The queue's `--replace invalid` cannot verify a `.fz` frame (`jobs.py` skips them).
- The raw archive sweep has not been run (Matthew, 2026-10-21: pilot before his leave, sweep after it).
