"""Put Astra raw frames that were transformed the wrong number of times back into the stored MaxIm orientation.

Calibration assumes Andor frames stored as 2048 x 2088 (NAXIS1 x NAXIS2), the orientation that
download/astra_transform.py produces with np.rot90(k=1), and SPIRIT frames mirrored left-right once
(ASTRAMIR='Y'). Audit of 2026-10-16:
  k1    Andor frame never rotated, ESO orientation (Callisto 2024: frames fetched late, four nights of calibrations)
  k3    Andor frame rotated twice, 180 deg from ESO (Ganymede 2025-11-15 and 11-25..29, Europa 2025-11-28)
  flip  SPIRIT frame mirrored twice, so back in ESO orientation with ASTRAMIR='Y' (Callisto 2025-12-22/23, ...);
        light frames show it as two 'WCS updated for horizontal mirroring' HISTORY cards
  mirror  SPIRIT frame never mirrored, in ESO orientation without ASTRAMIR (darks on Callisto 2025-04-24 and
        2026-03-05..09)
The operation is decided per frame outside this tool (the detector tells: hot pixels in darks and lights,
flat-field structure in flats and lights) and given to it in a CSV.

Report-only unless --apply. Input: a CSV with columns file,k where k is 1/k1, 3/k3, flip or mirror. A frame is
refused unless it is in the state its operation repairs: an Andor frame for k1/k3 must be in the native 2088 x 2048
shape; a frame for flip must be a 1024 x 1280 SPIRIT frame with ASTRAMIR='Y', no earlier repair card and, if it is
a light frame, at least two mirroring cards; a frame for mirror must be a 1024 x 1280 SPIRIT frame with no ASTRAMIR
and no mirroring card. So running the tool twice cannot repair a frame twice.

With --apply, per frame: the stored integers are turned with np.rot90(k) or flipped left-right (BITPIX, BZERO and
BSCALE kept); the celestial WCS is removed (it belonged to the old orientation; T7 with --force-platesolve solves
the frame again); for a turn ASTRAROT is set to 90 and IMAGEW/IMAGEH are swapped, for a flip ASTRAMIR stays 'Y', for a mirror
it is set to 'Y' as astra_transform would have;
HISTORY records the repair. The new file is written beside the old one under a hidden name and read back: the
stored integers and the physical values must equal the operation applied to the old ones. Then the old file is
hard-linked into --quarantine/<telescope>/<night>/ and the new file is renamed over the old name (atomic), and one
manifest row (old and new md5, operation, quarantine path) is appended. The quarantine must be on the same
filesystem as the frames.

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
SPIRIT_SHAPE = (1024, 1280)
REPAIR_CARD = 'orient repair'
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


def operation(value):
    """'1', 'k1', '3', 'k3' -> 1 or 3; 'flip', 'mirror' -> themselves; anything else -> None."""
    v = str(value).strip().lower()
    if v in ('flip', 'mirror'):
        return v
    v = v.lstrip('k')
    return int(v) if v in ('1', '3') else None


def apply_op(data, op):
    return np.flip(data, axis=1) if op in ('flip', 'mirror') else np.rot90(data, op)


def check(path, op):
    """None if the frame can be repaired by op, else the reason it is refused."""
    if op not in (1, 3, 'flip', 'mirror'):
        return f'refused: operation {op!r} (only k1, k3, flip and mirror are repairs)'
    h = fits.getheader(path)
    shape = (h.get('NAXIS1'), h.get('NAXIS2'))
    history = [str(x) for x in h.get('HISTORY', [])]
    if any(REPAIR_CARD in x for x in history):
        return 'refused: already repaired'
    if op == 'flip':
        if shape != SPIRIT_SHAPE:
            return f'refused: NAXIS {shape[0]}x{shape[1]} is not a SPIRIT frame'
        if str(h.get('ASTRAMIR', '')).strip() != 'Y':
            return 'refused: no ASTRAMIR=Y'
        cards = sum('horizontal mirroring' in x for x in history)
        if 'LIGHT' in str(h.get('IMAGETYP', '')).upper() and cards < 2:
            return f'refused: light frame with {cards} mirroring card(s), not mirrored twice'
        return None
    if op == 'mirror':
        if shape != SPIRIT_SHAPE:
            return f'refused: NAXIS {shape[0]}x{shape[1]} is not a SPIRIT frame'
        if 'ASTRAMIR' in h:
            return f"refused: ASTRAMIR={h['ASTRAMIR']!r}, already mirrored"
        cards = sum('horizontal mirroring' in x for x in history)
        if cards:
            return f'refused: {cards} mirroring card(s) without ASTRAMIR'
        return None
    if shape != NATIVE_SHAPE:
        return f'refused: NAXIS {shape[0]}x{shape[1]} is not the native shape'
    return None


def repaired_header(hdr, op, today):
    hdr = hdr.copy()
    for key in [c for c in list(hdr.keys()) if WCS_KEYS.match(c)]:
        del hdr[key]
    if op == 'flip':
        hdr.add_history(f'{today} {REPAIR_CARD}: fliplr; was mirrored twice')
    elif op == 'mirror':
        hdr['ASTRAMIR'] = ('Y', 'Y axis flip applied to astra images')
        hdr.add_history(f'{today} {REPAIR_CARD}: fliplr; was never mirrored')
    else:
        if 'IMAGEW' in hdr and 'IMAGEH' in hdr:
            hdr['IMAGEW'], hdr['IMAGEH'] = hdr['IMAGEH'], hdr['IMAGEW']
        hdr['ASTRAROT'] = (90, 'Anticlockwise rotation applied to astra images')
        was = 'never rotated' if op == 1 else 'rotated twice'
        hdr.add_history(f'{today} {REPAIR_CARD}: rot90 k={op}; was {was}')
    hdr.add_history(f'{today} {REPAIR_CARD}: celestial WCS removed; re-solve')
    return hdr


def repair(path, op, quarantine, today=None):
    today = today or time.strftime('%Y-%m-%d')
    with fits.open(path, do_not_scale_image_data=True, memmap=False) as hl:
        hdr = hl[0].header.copy()
        old = np.array(hl[0].data)
    new = np.ascontiguousarray(apply_op(old, op))
    bzero, bscale, bitpix = hdr.get('BZERO'), hdr.get('BSCALE'), hdr.get('BITPIX')
    shape = SPIRIT_SHAPE if op in ('flip', 'mirror') else STORED_SHAPE
    hdu = fits.PrimaryHDU(data=new, header=repaired_header(hdr, op, today), do_not_scale_image_data=True)
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
                     and h2.get('BSCALE') == bscale and (h2['NAXIS1'], h2['NAXIS2']) == shape)
    physical_ok = np.array_equal(fits.getdata(tmp).astype(np.float64),
                                 apply_op(fits.getdata(path).astype(np.float64), op))
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
    return dict(time=time.strftime('%Y-%m-%dT%H:%M:%S'), telescope=tel, night=night, file=name, k=op,
                old_md5=old_md5, new_md5=md5(path), bytes=os.path.getsize(path), quarantine=qpath)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('frames', help='CSV with columns file,k (k1, k3, flip or mirror)')
    ap.add_argument('--apply', action='store_true', help='repair the frames (default: report only)')
    ap.add_argument('--quarantine', help='directory for the originals, on the same filesystem as the frames')
    ap.add_argument('--manifest', help='CSV appended with one row per repaired frame')
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
        path, op = r['file'], operation(r['k'])
        reason = check(path, op)
        if reason:
            print(path, reason)
            refused += 1
            continue
        if not a.apply:
            print(f"would {op if op in ('flip', 'mirror') else f'turn by np.rot90(k={op})'} {path}")
            done += 1
            continue
        try:
            writer.writerow(repair(path, op, a.quarantine, today))
            mf.flush()
            done += 1
        except Exception as e:
            print(path, 'ERROR', e)
            refused += 1
    print(f"{'repaired' if a.apply else 'would repair'} {done}, refused or failed {refused}")
    return 1 if refused else 0


if __name__ == '__main__':
    sys.exit(main())
