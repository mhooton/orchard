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
            Fix: transform CRPIX and the matrix as astra_transform now does.  unmir is the same for a mirror.
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
HISTORY line, so a second run finds nothing to do.

Usage:  python repair_astra_wcs.py NIGHT_DIR [NIGHT_DIR ...] [--apply] [--csv report.csv]
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np
from astropy.io import fits

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from download.astra_transform import (MIRROR_X, ROT90_ACW, drop_singular_wcs, linear_matrix,  # noqa: E402
                                      plate_pa, set_linear_matrix, transform_wcs)

EXTS = ('*.fits', '*.fts', '*.fit')


def history(header):
    return ' '.join(str(h) for h in header.get('HISTORY', []))


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
    elif (sky and m is not None and not transformed_by_history and convention in ('CD', 'PC') and 'CDELT1' not in header
          and 'PA' in header and abs((float(header['PA']) - plate_pa(m) + 180) % 360 - 180) < 1.0):
        # The image was rotated or mirrored by an earlier downloader that left the whole WCS (CRPIX too) as it came
        # from ESO: no transform HISTORY, and the transfer's PA still agrees with the matrix. Seen on Callisto 2024.
        if rot:
            steps = int(round(float(rot) / 90.0)) % 4
            fixed = m @ np.linalg.matrix_power(np.array(ROT90_ACW, dtype=float), steps)
            found.append(('unrot', '{} x 90 deg'.format(steps)))
        else:
            fixed = m @ np.array(MIRROR_X, dtype=float)
            found.append(('unmir', 'mirror'))
    elif (sky and m is not None and transformed_by_history and 'PA' in header
          and abs((float(header['PA']) - plate_pa(m) + 180) % 360 - 180) < 1.0):
        # Transform HISTORY is there, yet the transfer's PA still agrees with the matrix, so the matrix is still the
        # native one; with more than one transform recorded (Ganymede 2025-07-28 has two), CRPIX cannot be trusted
        # either. Reported only: re-solve these frames instead.
        return [('ambiguous', 'matrix still native, CRPIX history unclear')]
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


def frames(night_dir):
    out = []
    for pattern in EXTS:
        out += glob.glob(os.path.join(night_dir, pattern))
    return sorted(out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    p.add_argument('night_dirs', nargs='+')
    p.add_argument('--apply', action='store_true', help='write the fixes (default: report only)')
    p.add_argument('--csv', help='write one row per affected frame here')
    args = p.parse_args(argv)
    rows, counts = [], {}
    for nd in args.night_dirs:
        for path in frames(nd):
            try:
                with fits.open(path, mode='update' if args.apply else 'readonly', memmap=False) as h:
                    found = diagnose(h[0].header)
                    if not found:
                        continue
                    for k, _ in found:
                        counts[k] = counts.get(k, 0) + 1
                    rows.append((path, ';'.join('{}:{}'.format(k, d) for k, d in found)))
                    if args.apply and not any(k in ('image_not_rotated', 'ambiguous') for k, _ in found):
                        repair(h[0].header, found)
            except Exception as e:
                rows.append((path, 'error:{}'.format(e)))
                counts['error'] = counts.get('error', 0) + 1
    if args.csv:
        with open(args.csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['path', 'found'])
            w.writerows(rows)
    print('{} frames affected ({}); {}'.format(len(rows), ', '.join('{} {}'.format(v, k) for k, v in sorted(counts.items())) or 'none',
                                              'repaired' if args.apply else 'nothing written (report only)'))


if __name__ == '__main__':
    main()
