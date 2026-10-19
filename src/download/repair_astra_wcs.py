"""
repair_astra_wcs.py - correct the WCS that astra_transform left in Astra raw frames before ASTRAWCS existed.

Nothing is written unless --apply is given. Only raw light frames are considered, and only these signatures:

  rot_pc    ASTRAROT set, a PCi_j matrix with no CDELTi, the transform's 'WCS updated for ... rotation' HISTORY,
            no ASTRAWCS. The image was rotated and CRPIX moved, but the matrix is still the one Astra wrote for
            the native frame, so the WCS is 90 degrees out about CRPIX.  Fix: PC <- PC @ ROT90_ACW per 90 degrees.
  mir_pc    The same for SPIRIT: ASTRAMIR = 'Y', PC without CDELT, the 'mirroring' HISTORY.  Fix: PC <- PC @ MIRROR_X.
  singular  A sky WCS whose matrix is singular (SPIRIT frames arrive with an all-zero CD block and PA = 180.0).
            Fix: remove the WCS and PA, so astrometry.py solves the frame.
  unrot     ASTRAROT set but no transform HISTORY, and the transfer's PA still agrees with the matrix: an earlier
            downloader rotated the image and left the whole WCS, CRPIX included, as ESO delivered it (Callisto 2024).
            Fix: transform CRPIX and the matrix as astra_transform now does.
  unmir     The same for a mirrored SPIRIT frame, detected by parity (a mirror leaves PA unchanged): no transform
            HISTORY and a matrix with the ESO parity (-1; every stored SPIRIT frame solved independently has +1).
  stale_pa  Any transformed frame whose PA describes the native orientation rather than the stored one.
  image_not_rotated  ASTRAROT on an Andor image still 2088 x 2048: the WCS was turned and the image was not
            (Ganymede, November 2025). Reported only.
  ambiguous Transform HISTORY present but the matrix still agrees with the transfer's PA, so it was never turned
            (or turned and turned back); CRPIX may have been moved more than once (Ganymede 2025-07-28). Reported
            only; these frames need a fresh plate solve.
            Fix: recompute PA from the matrix as it is used (wcslib precedence: PC over CD).

A PC matrix written with CDELTi = 1 is orchard's own plate solution (astropy writes it that way), made on the
stored image, and is left alone; so is a CD matrix (astrometry.net's, rotated by the old CD formula, which is
within ~1e-6 deg/px of right for these nearly skew-free matrices). Repaired frames get ASTRAWCS = 2 and a
HISTORY line, so a second run finds nothing to do, and the edit record of every in-place edit (ORCHMOD and a card
'orchard <version> repair_astra_wcs: <signatures>', docs/versioning.md).

Writing touches header bytes only: the data are never decoded, the original header of every changed frame is saved
first (one tar per night), and file mode and timestamps are kept. --restore puts the saved headers back. CHECKSUM
is recomputed from the new header and the frame's DATASUM, which is trusted when the old CHECKSUM agrees with it
and otherwise computed from the data bytes (absent, or stale since astra_transform rotated the image), so every
repaired frame verifies.

Usage:  python repair_astra_wcs.py NIGHT_DIR [NIGHT_DIR ...] [--apply --backup DIR] [--csv report.csv]
        python repair_astra_wcs.py --restore DIR/<telescope>_images_<night>.tar
"""
import argparse
import csv
import glob
import hashlib
import io
import os
import shutil
import sys
import tarfile

import numpy as np
from astropy.io import fits

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from download.astra_transform import (MIRROR_X, ROT90_ACW, drop_singular_wcs, linear_matrix,  # noqa: E402
                                      plate_pa, set_linear_matrix, transform_wcs)
from download.raw_provenance import claimed_datasum, data_sum, set_checksum  # noqa: E402
from utils import orchard_version  # noqa: E402

EXTS = ('*.fits', '*.fts', '*.fit')


def history(header):
    return ' '.join(str(h) for h in header.get('HISTORY', []))


def pa_agrees(header, m, tol=1.0):
    """True if the transfer's PA keyword still describes matrix m (within tol degrees)."""
    if 'PA' not in header:
        return False
    return abs((float(header['PA']) - plate_pa(m) + 180) % 360 - 180) < tol


def diagnose(header):
    """List of (signature, detail) for one header. Pure: does not modify the header."""
    if str(header.get('IMAGETYP', '')) != 'Light Frame' or 'ASTRAWCS' in header:
        return []
    rot = header.get('ASTRAROT')
    mirrored = str(header.get('ASTRAMIR', '')).strip() == 'Y'
    if not rot and not mirrored:
        return []
    found = []
    m, convention = linear_matrix(header)
    sky = str(header.get('CTYPE1', '')).startswith('RA')
    if rot and (header.get('NAXIS1'), header.get('NAXIS2')) == (2088, 2048):
        # ASTRAROT set on an Andor image that is still in the native shape (the 'rotate_wcs_only' mode, seen on
        # Ganymede in November 2025): the WCS was turned and the image was not. Reported, never fixed blind.
        return [('image_not_rotated', 'native 2088 x 2048 with ASTRAROT')]
    if sky and m is not None and abs(np.linalg.det(m)) == 0:
        return [('singular', 'matrix is singular')]
    pc_without_cdelt = convention == 'PC' and 'CDELT1' not in header and 'CDELT2' not in header
    hist = history(header)
    transformed_by_history = 'anticlockwise rotation' in hist or 'horizontal mirroring' in hist
    fixed = m
    if sky and pc_without_cdelt and rot and 'anticlockwise rotation' in hist:
        steps = int(round(float(rot) / 90.0)) % 4
        fixed = m @ np.linalg.matrix_power(np.array(ROT90_ACW, dtype=float), steps)
        found.append(('rot_pc', '{} x 90 deg'.format(steps)))
    elif sky and pc_without_cdelt and mirrored and 'horizontal mirroring' in hist:
        fixed = m @ np.array(MIRROR_X, dtype=float)
        found.append(('mir_pc', 'mirror'))
    elif (rot and sky and m is not None and not transformed_by_history and convention in ('CD', 'PC')
          and 'CDELT1' not in header and pa_agrees(header, m)):
        # The image was rotated by an earlier downloader that left the whole WCS (CRPIX too) as it came from ESO:
        # no transform HISTORY, and the transfer's PA still agrees with the matrix. Seen on Callisto in 2024.
        steps = int(round(float(rot) / 90.0)) % 4
        fixed = m @ np.linalg.matrix_power(np.array(ROT90_ACW, dtype=float), steps)
        found.append(('unrot', '{} x 90 deg'.format(steps)))
    elif rot and sky and m is not None and transformed_by_history and pa_agrees(header, m):
        # Transform HISTORY is there, yet the transfer's PA still agrees with the matrix, so the matrix is still the
        # native one; with more than one transform recorded (Ganymede 2025-07-28 has two), CRPIX cannot be trusted
        # either. Reported only: re-solve these frames instead.
        return [('ambiguous', 'matrix still native, CRPIX history unclear')]
    elif mirrored and sky and m is not None and np.linalg.det(m) < 0:
        # A mirror leaves PA (the angle of -Y) unchanged, so for SPIRIT the test is parity: every SPIRIT frame we
        # store, plate-solved independently, has parity +1, and the ESO original has -1.
        if transformed_by_history or 'CDELT1' in header:
            # either the transform ran and the matrix is still native, or the matrix is orchard's own solution (PC
            # with CDELT = 1, made on the stored image) and came out mirrored; neither is a header fix
            return [('ambiguous', 'mirrored image with the native parity')]
        fixed = m @ np.array(MIRROR_X, dtype=float)
        found.append(('unmir', 'mirror'))
    if sky and fixed is not None and 'PA' in header:
        pa = plate_pa(fixed)
        if abs((float(header['PA']) - pa + 180) % 360 - 180) > 0.01:
            found.append(('stale_pa', '{:.3f} -> {:.3f}'.format(float(header['PA']), pa)))
    return found


def repair(header, found):
    """Apply the fixes listed by diagnose() to header in place."""
    kinds = {k for k, _ in found}
    if kinds & {'image_not_rotated', 'ambiguous'}:
        return
    if 'singular' in kinds:
        drop_singular_wcs(header)
    elif 'unrot' in kinds:
        # The whole WCS is still the native one: move CRPIX as well, starting from the native dimensions
        steps = int(round(float(header['ASTRAROT']) / 90.0)) % 4
        nx, ny = header['NAXIS1'], header['NAXIS2']
        if steps % 2:
            nx, ny = ny, nx
        for _ in range(steps):
            transform_wcs(header, ROT90_ACW, lambda c1, c2, nx=nx: (c2, nx - c1 + 1))
            nx, ny = ny, nx
    elif 'unmir' in kinds:
        nx = header['NAXIS1']
        transform_wcs(header, MIRROR_X, lambda c1, c2: (nx - c1 + 1, c2))
    else:
        m, convention = linear_matrix(header)
        if 'rot_pc' in kinds:
            steps = int(round(float(header['ASTRAROT']) / 90.0)) % 4
            m = m @ np.linalg.matrix_power(np.array(ROT90_ACW, dtype=float), steps)
            set_linear_matrix(header, m, convention)
        elif 'mir_pc' in kinds:
            m = m @ np.array(MIRROR_X, dtype=float)
            set_linear_matrix(header, m, convention)
        if 'PA' in header and m is not None:
            header['PA'] = (plate_pa(m), '[deg, 0-360 CCW] Position angle of plate')
    header['ASTRAWCS'] = (2, 'WCS matrix transformed with the image')
    header.add_history('repair_astra_wcs: ' + ', '.join(sorted(kinds)))
    orchard_version.record_edit(header, 'repair_astra_wcs', ', '.join(sorted(kinds)))


def frames(night_dir):
    out = []
    for pattern in EXTS:
        out += glob.glob(os.path.join(night_dir, pattern))
    return sorted(out)


# Writing. These are the archive's raw frames, so the data are never decoded: the header is read and written as
# bytes, the data bytes stay where they are, and the original header goes to a tar file first so the repair can be
# undone with --restore. Decoding scaled 16-bit data through astropy's update mode can rewrite it as float.
BLOCK = 2880


def read_header_bytes(path):
    """The primary header exactly as stored (to the END card, padded to whole blocks), and its parse."""
    chunks = []
    with open(path, 'rb') as f:
        while True:
            block = f.read(BLOCK)
            if len(block) < BLOCK:
                raise ValueError('{}: no END card before the end of the file'.format(path))
            chunks.append(block)
            if any(block[i:i + 8] == b'END     ' for i in range(0, BLOCK, 80)):
                raw = b''.join(chunks)
                return raw, fits.Header.fromstring(raw.decode('ascii'))


def write_header_bytes(path, new_raw, old_raw):
    """
    Replace the header of path. In place when the new header has as many blocks as the old one, which is the
    usual case; otherwise via a temporary file holding the new header and the untouched data bytes, swapped in
    with os.replace. Keeps the file's mode and timestamps.
    """
    st = os.stat(path)
    if len(new_raw) == len(old_raw):
        with open(path, 'r+b') as f:
            f.seek(0)
            f.write(new_raw)
            f.flush()
            os.fsync(f.fileno())
    else:
        tmp = path + '.repair-tmp'
        with open(path, 'rb') as src, open(tmp, 'wb') as dst:
            src.seek(len(old_raw))
            dst.write(new_raw)
            shutil.copyfileobj(src, dst, 1 << 20)
            dst.flush()
            os.fsync(dst.fileno())
        if os.path.getsize(tmp) != st.st_size - len(old_raw) + len(new_raw):
            os.remove(tmp)
            raise IOError('{}: size check failed, left untouched'.format(path))
        os.chmod(tmp, st.st_mode & 0o7777)
        os.replace(tmp, path)
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns))
    if read_header_bytes(path)[0] != new_raw:
        raise IOError('{}: header read back differs from what was written'.format(path))


def data_digest(path, header_len):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        f.seek(header_len)
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def apply_to_file(path, backup):
    """Diagnose path and repair it on disk. Returns the list from diagnose(); backup is an open tarfile."""
    old_raw, header = read_header_bytes(path)
    found = diagnose(header)
    if not found or any(k in ('image_not_rotated', 'ambiguous') for k, _ in found):
        return found
    datasum = claimed_datasum(old_raw, header)
    if datasum is None:
        datasum = data_sum(path, len(old_raw), header)
    repair(header, found)
    set_checksum(header, datasum)
    new_raw = header.tostring(padding=True).encode('ascii')
    info = tarfile.TarInfo(os.path.abspath(path).lstrip('/') + '.hdr')
    info.size = len(old_raw)
    info.mtime = os.stat(path).st_mtime
    backup.addfile(info, io.BytesIO(old_raw))
    backup.fileobj.flush()
    if len(new_raw) == len(old_raw):
        # written in place over the header blocks only; the size check and the read-back in write_header_bytes
        # cover it, and reading every frame's data twice would cost ~14 TB over the archive
        size = os.path.getsize(path)
        write_header_bytes(path, new_raw, old_raw)
        if os.path.getsize(path) != size:
            raise IOError('{}: file size changed'.format(path))
    else:
        before = data_digest(path, len(old_raw))
        write_header_bytes(path, new_raw, old_raw)
        if data_digest(path, len(new_raw)) != before:
            raise IOError('{}: data bytes changed'.format(path))
    return found


def restore(tar_path):
    """Put back every header saved in tar_path. Returns the number of files restored."""
    n = 0
    with tarfile.open(tar_path) as t:
        for m in t.getmembers():
            path = '/' + m.name[:-len('.hdr')]
            old_raw = t.extractfile(m).read()
            current_raw, _ = read_header_bytes(path)
            write_header_bytes(path, old_raw, current_raw)
            n += 1
    return n


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    p.add_argument('night_dirs', nargs='*')
    p.add_argument('--apply', action='store_true', help='write the fixes (default: report only)')
    p.add_argument('--backup', help='directory for the per-night tar of original headers (required with --apply)')
    p.add_argument('--restore', help='put back the headers saved in this tar file, and do nothing else')
    p.add_argument('--csv', help='write one row per affected frame here')
    args = p.parse_args(argv)
    if args.restore:
        print('{} headers restored from {}'.format(restore(args.restore), args.restore))
        return
    if args.apply and not args.backup:
        p.error('--apply needs --backup')
    rows, counts = [], {}
    for nd in args.night_dirs:
        backup = None
        if args.apply:
            os.makedirs(args.backup, exist_ok=True)
            name = '_'.join(os.path.abspath(nd).rstrip('/').split('/')[-3:]) + '.tar'
            backup = tarfile.open(os.path.join(args.backup, name), 'a')
        try:
            for path in frames(nd):
                try:
                    if args.apply:
                        found = apply_to_file(path, backup)
                    else:
                        found = diagnose(read_header_bytes(path)[1])
                except Exception as e:
                    rows.append((path, 'error:{}'.format(e)))
                    counts['error'] = counts.get('error', 0) + 1
                    continue
                if not found:
                    continue
                for k, _ in found:
                    counts[k] = counts.get(k, 0) + 1
                rows.append((path, ';'.join('{}:{}'.format(k, d) for k, d in found)))
        finally:
            if backup is not None:
                backup.close()
    if args.csv:
        with open(args.csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['path', 'found'])
            w.writerows(rows)
    print('{} frames affected ({}); {}'.format(len(rows), ', '.join('{} {}'.format(v, k) for k, v in sorted(counts.items())) or 'none',
                                              'repaired' if args.apply else 'nothing written (report only)'))


if __name__ == '__main__':
    main()
