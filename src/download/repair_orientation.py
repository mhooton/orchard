"""Turn Andor raw frames that are in the native (ESO) shape into the stored MaxIm orientation.

Calibration assumes every Andor frame is stored as 2048 x 2088 (NAXIS1 x NAXIS2), the orientation that
download/astra_transform.py produces with np.rot90(k=1). Two kinds of frame ended up in the native 2088 x 2048
shape instead (audit of 2026-10-16):
  k=1  never rotated, ESO orientation (Callisto 2024: frames fetched late, and four nights of calibrations)
  k=3  rotated twice, 180 deg from ESO (Ganymede 2025-11-15 and 11-25..29, turned again on 2025-12-09)
The shape cannot tell them apart; the detector can (hot pixels in darks and lights, flat-field structure in flats
and lights), so k is decided per frame outside this tool and given to it in a CSV.

Report-only unless --apply. Input: a CSV with columns file,k (k may be written k1/k3). A frame that is not in the
native shape is refused, so running the tool twice cannot turn a frame twice.

With --apply, per frame: the stored integers are turned with np.rot90(k) (BITPIX, BZERO and BSCALE kept);
ASTRAROT is set to 90; the celestial WCS is removed (it belonged to the old orientation; T7 with
--force-platesolve solves the frame again); IMAGEW/IMAGEH are swapped; HISTORY records the repair. The new file is
written beside the old one under a hidden name and read back: the stored integers must equal rot90(old, k) and the
physical values rot90(old physical values, k). Then the old file is hard-linked into
--quarantine/<telescope>/<night>/ and the new file is renamed over the old name (atomic), and one manifest row
(old and new md5, k, quarantine path) is appended. The quarantine must be on the same filesystem as the frames.

    python download/repair_orientation.py frames.csv
    python download/repair_orientation.py frames.csv --apply --quarantine DIR --manifest manifest.csv
"""
import argparse
import csv
import hashlib
import os
import re
import sys
import time

import numpy as np
from astropy.io import fits

NATIVE_SHAPE = (2088, 2048)    # (NAXIS1, NAXIS2) as ESO delivers an Andor frame
STORED_SHAPE = (2048, 2088)
WCS_KEYS = re.compile(r'^(WCSAXES|CTYPE\d|CUNIT\d|CRVAL\d|CRPIX\d|CDELT\d|CROTA\d|CD\d_\d|PC\d_\d|PV\d_\d+|'
                      r'LONPOLE|LATPOLE|RADESYS|EQUINOX|A_ORDER|B_ORDER|AP_ORDER|BP_ORDER|A_\d_\d|B_\d_\d|'
                      r'AP_\d_\d|BP_\d_\d|A_DMAX|B_DMAX)$')
MANIFEST_FIELDS = ['time', 'telescope', 'night', 'file', 'k', 'old_md5', 'new_md5', 'bytes', 'quarantine']


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def telescope_night(path):
    m = re.search(r'/Observations/([^/]+)/images/([^/]+)/', path)
    return (m.group(1), m.group(2)) if m else ('unknown', 'unknown')


def check(path, k):
    """None if the frame can be turned by k, else the reason it is refused."""
    if k not in (1, 3):
        return f'refused: k={k} (only 1 and 3 turn a native frame into the stored shape)'
    h = fits.getheader(path)
    if (h.get('NAXIS1'), h.get('NAXIS2')) != NATIVE_SHAPE:
        return f"refused: NAXIS {h.get('NAXIS1')}x{h.get('NAXIS2')} is not the native shape"
    return None


def turned_header(hdr, k, today):
    hdr = hdr.copy()
    for key in [c for c in list(hdr.keys()) if WCS_KEYS.match(c)]:
        del hdr[key]
    if 'IMAGEW' in hdr and 'IMAGEH' in hdr:
        hdr['IMAGEW'], hdr['IMAGEH'] = hdr['IMAGEH'], hdr['IMAGEW']
    hdr['ASTRAROT'] = (90, 'Anticlockwise rotation applied to astra images')
    was = 'never rotated' if k == 1 else 'rotated twice'
    hdr.add_history(f'{today} orient repair: rot90 k={k}; was {was}')
    hdr.add_history(f'{today} orient repair: celestial WCS removed; re-solve')
    return hdr


def repair(path, k, quarantine, today=None):
    today = today or time.strftime('%Y-%m-%d')
    with fits.open(path, do_not_scale_image_data=True, memmap=False) as hl:
        hdr = hl[0].header.copy()
        old = np.array(hl[0].data)
    new = np.ascontiguousarray(np.rot90(old, k))
    bzero, bscale, bitpix = hdr.get('BZERO'), hdr.get('BSCALE'), hdr.get('BITPIX')
    hdu = fits.PrimaryHDU(data=new, header=turned_header(hdr, k, today), do_not_scale_image_data=True)
    # PrimaryHDU(do_not_scale_image_data=True) drops BZERO/BSCALE: put them back as they were
    if bzero is not None:
        hdu.header['BZERO'] = bzero
    if bscale is not None:
        hdu.header['BSCALE'] = bscale
    directory, name = os.path.split(path)
    tmp = os.path.join(directory, f'.{name}.orient_repair.tmp')
    hdu.writeto(tmp, overwrite=True, output_verify='silentfix')
    with fits.open(tmp, do_not_scale_image_data=True, memmap=False) as chk:
        h2 = chk[0].header
        stored_ok = (np.array_equal(chk[0].data, new) and h2.get('BITPIX') == bitpix and h2.get('BZERO') == bzero
                     and h2.get('BSCALE') == bscale and (h2['NAXIS1'], h2['NAXIS2']) == STORED_SHAPE)
    physical_ok = np.array_equal(fits.getdata(tmp).astype(np.float64),
                                 np.rot90(fits.getdata(path).astype(np.float64), k))
    if not (stored_ok and physical_ok):
        os.unlink(tmp)
        raise RuntimeError(f'verification failed (stored values {stored_ok}, physical values {physical_ok})')
    tel, night = telescope_night(path)
    qdir = os.path.join(quarantine, tel, night)
    os.makedirs(qdir, exist_ok=True)
    qpath = os.path.join(qdir, name)
    old_md5 = md5(path)
    os.link(path, qpath)
    os.replace(tmp, path)
    return dict(time=time.strftime('%Y-%m-%dT%H:%M:%S'), telescope=tel, night=night, file=name, k=k,
                old_md5=old_md5, new_md5=md5(path), bytes=os.path.getsize(path), quarantine=qpath)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('frames', help='CSV with columns file,k')
    ap.add_argument('--apply', action='store_true', help='turn the frames (default: report only)')
    ap.add_argument('--quarantine', help='directory for the originals, on the same filesystem as the frames')
    ap.add_argument('--manifest', help='CSV appended with one row per turned frame')
    a = ap.parse_args(argv)
    if a.apply and not (a.quarantine and a.manifest):
        ap.error('--apply needs --quarantine and --manifest')
    rows = list(csv.DictReader(open(a.frames)))
    today = time.strftime('%Y-%m-%d')
    writer = None
    if a.apply:
        new = not os.path.exists(a.manifest)
        mf = open(a.manifest, 'a', newline='')
        writer = csv.DictWriter(mf, fieldnames=MANIFEST_FIELDS)
        if new:
            writer.writeheader()
    done = refused = 0
    for r in rows:
        path, k = r['file'], int(str(r['k']).lstrip('k'))
        reason = check(path, k)
        if reason:
            print(path, reason)
            refused += 1
            continue
        if not a.apply:
            print(f'would turn {path} by np.rot90(k={k})')
            done += 1
            continue
        try:
            writer.writerow(repair(path, k, a.quarantine, today))
            mf.flush()
            done += 1
        except Exception as e:
            print(path, 'ERROR', e)
            refused += 1
    print(f"{'turned' if a.apply else 'would turn'} {done}, refused or failed {refused}")
    return 1 if refused else 0


if __name__ == '__main__':
    sys.exit(main())
